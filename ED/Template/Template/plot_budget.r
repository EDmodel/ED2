#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"     # Current directory.
srcdir         = "thisrscpath"  # Source  directory.
outroot        = "thisoutroot"  # Output  directory.
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

outform        = thisoutform           # Formats for output file.  Supported formats are:
                                 #   - "X11" - for printing on screen
                                 #   - "eps" - for postscript printing
                                 #   - "png" - for PNG printing
                                 #   - "pdf" - for PDF printing

cex.main       = 0.8             # Scale coefficient for the title

byeold         = TRUE           # Remove old files of the given format?

depth          = 96             # PNG resolution, in pixels per inch
paper          = "square"       # Paper size, to define the plot shape
ptsz           = 22             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.05           # inset distance between legend and edge of plot region.

scalleg        = 0.32           # Increase in y scale to fit the legend.
ncolshov       = 200            # Target number of colours for Hovmoller diagrams.
plotgrid       = FALSE           # Should I include a grid in plots?
ibackground    = mybackground   # Background settings (check load_everything.r)
f.leg          = 1/6            # Fraction of plotting area for legend.


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#


#----- Time series plots. -----------------------------------------------------------------#
n           = 0
budget      = list()
n           = n + 1
budget[[n]] = list( vnam   = c("co2.nep.rel","co2.eddy.flux.rel"
                              ,"co2.dstorage.rel","co2.residual.rel")
                  , desc   = c("NEP","Eddy flux","Delta (Storage)","Residual")
                  , colour = c("#A3CC52","#2996CC","#990F0F","grey30")
                  , lty    = c("longdash","twodash","dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "co2.rel"
                  , theme  = "Accumulated Carbon dioxide (relative)"
                  , unit   = "pc"
                  , ylim   = NA
                  , aggr   = FALSE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("ene.precip.rel","ene.netrad.rel","ene.prss.eff.rel"
                              ,"ene.eddy.flux.rel","ene.runoff.rel","ene.dstorage.rel"
                              ,"ene.residual.rel")
                  , desc   = c("Rainfall","Net Radiation","Pressure effect","Eddy flux"
                              ,"Total runoff","Delta (Storage)","Residual")
                  , colour = c("#3B24B3","#E65C17","#306614","#2996CC"
                              ,"#A3CC52","#990F0F","grey30")
                  , lty    = c("longdash","dotdash","dashed","twodash"
                              ,"longdash"  ,"dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "ene.rel"
                  , theme  = "Accumulated enthalpy (relative)"
                  , unit   = "pc"
                  , ylim   = c(-1.5,1.5)
                  , aggr   = FALSE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("h2o.precip.rel","h2o.eddy.flux.rel","h2o.runoff.rel"
                              ,"h2o.dstorage.rel","h2o.residual.rel")
                  , desc   = c("Rainfall","Eddy flux","Total Runoff"
                              ,"Delta (Storage)","Residual")
                  , colour = c("#3B24B3","#2996CC","#A3CC52","#990F0F","grey30")
                  , lty    = c("longdash","twodash","longdash","dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "h2o.rel"
                  , theme  = "Accumulated water (relative)"
                  , unit   = "pc"
                  , ylim   = c(-1.5,1.5)
                  , aggr   = FALSE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("co2.nep.cum","co2.eddy.flux.cum"
                              ,"co2.dstorage.cum","co2.residual.cum")
                  , desc   = c("NEP","Eddy flux","Delta (Storage)","Residual")
                  , colour = c("#A3CC52","#2996CC","#990F0F","grey30")
                  , lty    = c("longdash","twodash","dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "co2.cum"
                  , theme  = "Accumulated Carbon dioxide"
                  , unit   = "umolcom2"
                  , ylim   = NA
                  , aggr   = FALSE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("ene.precip.cum","ene.netrad.cum","ene.prss.eff.cum"
                              ,"ene.eddy.flux.cum","ene.runoff.cum","ene.dstorage.cum"
                              ,"ene.residual.cum")
                  , desc   = c("Rainfall","Net Radiation","Pressure effect","Eddy flux"
                              ,"Total runoff","Delta (Storage)","Residual")
                  , colour = c("#3B24B3","#E65C17","#306614","#2996CC"
                              ,"#A3CC52","#990F0F","grey30")
                  , lty    = c("longdash","dotdash","dashed","twodash"
                              ,"longdash"  ,"dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "ene.cum"
                  , theme  = "Accumulated enthalpy"
                  , unit   = "jom2"
                  , ylim   = NA
                  , aggr   = FALSE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("h2o.precip.cum","h2o.eddy.flux.cum","h2o.runoff.cum"
                              ,"h2o.dstorage.cum","h2o.residual.cum")
                  , desc   = c("Rainfall","Eddy flux","Total Runoff"
                              ,"Delta (Storage)","Residual")
                  , colour = c("#3B24B3","#2996CC","#A3CC52","#990F0F","grey30")
                  , lty    = c("longdash","twodash","longdash","dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "h2o.cum"
                  , theme  = "Accumulated water"
                  , unit   = "kgwom2"
                  , ylim   = NA
                  , aggr   = FALSE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("co2.nep","co2.eddy.flux","co2.dstorage","co2.residual")
                  , desc   = c("NEP","Eddy flux","Delta (Storage)","Residual")
                  , colour = c("#A3CC52" ,"#2996CC","#990F0F","grey30")
                  , lty    = c("longdash","twodash","dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0)
                  , range  = c(TRUE,TRUE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "carbflux"
                  , theme  = "Carbon dioxide budget"
                  , unit   = "umolom2"
                  , ylim   = NA
                  , aggr   = TRUE
                  , plt    = TRUE
                  )#end list
n           = n + 1
budget[[n]] = list( vnam   = c("ene.precip","ene.netrad","ene.prss.eff","ene.eddy.flux"
                              ,"ene.runoff","ene.dstorage","ene.residual")
                  , desc   = c("Rainfall","Net Radiation","Pressure effect","Eddy flux"
                              ,"Total runoff","Delta (Storage)","Residual")
                  , colour = c("#3B24B3","#E65C17","#306614","#2996CC"
                              ,"#A3CC52","#990F0F","grey30")
                  , lty    = c("longdash","dotdash","dashed","twodash"
                              ,"longdash"  ,"dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0,3.0,3.0,3.0)
                  , range  = c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "eneflux"
                  , theme  = "Enthalpy budget"
                  , unit   = "wom2"
                  , ylim   = NA
                  , aggr   = TRUE
                  , plt    = TRUE)
n           = n + 1
budget[[n]] = list( vnam   = c("h2o.precip","h2o.eddy.flux","h2o.runoff"
                              ,"h2o.dstorage","h2o.residual")
                  , desc   = c("Rainfall","Eddy flux","Total Runoff"
                              ,"Delta (Storage)","Residual")
                  , colour = c("#3B24B3","#2996CC","#A3CC52","#990F0F","grey30")
                  , lty    = c("longdash","twodash","longdash","dotdash","solid")
                  , lwd    = c(3.0,3.0,3.0,3.0,3.0)
                  , range  = c(FALSE,TRUE,FALSE,FALSE,TRUE)
                  , type   = ptype
                  , plog   = ""
                  , prefix = "h2oflux"
                  , theme  = "Water budget"
                  , unit   = "kgwom2oday"
                  , ylim   = NA
                  , aggr   = TRUE
                  , plt    = TRUE)
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
ex.size = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
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
   inpref     = file.path(here,place)
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
   mypatches = length(grep("budget_state_patch_",filelist))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Patch loop.                                                                        #
   #---------------------------------------------------------------------------------------#
   for (ipa in sequence(mypatches)){
      #----- Find the character version of the patch number. ------------------------------#
      cipa = sprintf("%4.4i",ipa)

      cat0("    - Patch # ",ipa,".")
      #----- Define the output directory. -------------------------------------------------#
      patchdir  = file.path(outdir,paste("patch_",cipa,sep=""))
      if (! file.exists(patchdir)) dir.create(patchdir)

      #----- Define the input file name. --------------------------------------------------#
      inputfile = file.path(inpref,paste("budget_state_patch_",cipa,".txt",sep=""))
      cat0("      * Open file:",inputfile)

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

      #----- Reduce the size of the file to be the period of interest only. ---------------#
      cat0("      * Reduce data to the period of interest.")
      when   = chron( chron(dates=paste(cpatch$month,cpatch$day,cpatch$year,sep="/"))
                    + cpatch$time/day.sec, out.format=c(dates="m/d/y",times="h:m:s"))
      sel    = when >= whena & when <= whenz
      cpatch = cpatch[sel,]
      when   = when[sel]


      #------------------------------------------------------------------------------------#
      #     Reduce the number of variables, and standardise signs:                         #
      # > 0 : entering the system;                                                         #
      # < 0 : leaving the system.                                                          #
      #------------------------------------------------------------------------------------#
      cat0("      * Reduce the number of terms.")
      #----- CO2. -------------------------------------------------------------------------#
      cpatch$co2.eddy.flux = - cpatch$co2.loss2atm
      cpatch$co2.residual  =   cpatch$co2.residual + cpatch$co2.dens.eff
      #----- Enthalpy. --------------------------------------------------------------------#
      cpatch$ene.eddy.flux = - cpatch$ene.loss2atm
      cpatch$ene.runoff    = - cpatch$ene.runoff   - cpatch$ene.drainage
      cpatch$ene.residual  =   cpatch$ene.residual + cpatch$ene.dens.eff
      #----- Water. -----------------------------------------------------------------------#
      cpatch$h2o.eddy.flux = - cpatch$h2o.loss2atm
      cpatch$h2o.runoff    = - cpatch$h2o.runoff   - cpatch$h2o.drainage
      cpatch$h2o.residual  =   cpatch$h2o.residual + cpatch$h2o.dens.eff
      #------------------------------------------------------------------------------------#


      #----- Re-scale or re-define some variables. ----------------------------------------#
      cat0("      * Find the average storage.")
      co2.scale = mean(cpatch$co2.storage)
      ene.scale = mean(cpatch$ene.storage)
      h2o.scale = mean(cpatch$h2o.storage)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Re-scale or re-define some variables.                                         #
      #------------------------------------------------------------------------------------#
      cat0("      * Define the cumulative sum of all budget terms.")
      #----- Cumulative terms for all CO2 budget. -----------------------------------------#
      cpatch$co2.dstorage.cum  = cumsum(cpatch$co2.dstorage )
      cpatch$co2.nep.cum       = cumsum(cpatch$co2.nep      )
      cpatch$co2.eddy.flux.cum = cumsum(cpatch$co2.eddy.flux)
      cpatch$co2.residual.cum  = cumsum(cpatch$co2.residual )
      #----- Cumulative terms for all enthalpy budget. ------------------------------------#
      cpatch$ene.dstorage.cum  = cumsum(cpatch$ene.dstorage )
      cpatch$ene.precip.cum    = cumsum(cpatch$ene.precip   )
      cpatch$ene.netrad.cum    = cumsum(cpatch$ene.netrad   )
      cpatch$ene.prss.eff.cum  = cumsum(cpatch$ene.prss.eff )
      cpatch$ene.eddy.flux.cum = cumsum(cpatch$ene.eddy.flux)
      cpatch$ene.runoff.cum    = cumsum(cpatch$ene.runoff   )
      cpatch$ene.residual.cum  = cumsum(cpatch$ene.residual )
      #----- Cumulative terms for all water budget. ---------------------------------------#
      cpatch$h2o.dstorage.cum  = cumsum(cpatch$h2o.dstorage ) / day.sec
      cpatch$h2o.precip.cum    = cumsum(cpatch$h2o.precip   ) / day.sec
      cpatch$h2o.eddy.flux.cum = cumsum(cpatch$h2o.eddy.flux) / day.sec
      cpatch$h2o.runoff.cum    = cumsum(cpatch$h2o.runoff   ) / day.sec
      cpatch$h2o.residual.cum  = cumsum(cpatch$h2o.residual ) / day.sec
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Re-scale or re-define some variables.                                         #
      #------------------------------------------------------------------------------------#
      cat0("      * Define the relative cumulative sum of the terms.")
      #----- Cumulative terms for all CO2 budget. -----------------------------------------#
      cpatch$co2.dstorage.rel  = 100. * cpatch$co2.dstorage.cum  / co2.scale
      cpatch$co2.nep.rel       = 100. * cpatch$co2.nep.cum       / co2.scale
      cpatch$co2.eddy.flux.rel = 100. * cpatch$co2.eddy.flux.cum / co2.scale
      cpatch$co2.residual.rel  = 100. * cpatch$co2.residual.cum  / co2.scale
      #----- Cumulative terms for all enthalpy budget. ------------------------------------#
      cpatch$ene.dstorage.rel  = 100. * cpatch$ene.dstorage.cum  / ene.scale
      cpatch$ene.precip.rel    = 100. * cpatch$ene.precip.cum    / ene.scale
      cpatch$ene.netrad.rel    = 100. * cpatch$ene.netrad.cum    / ene.scale
      cpatch$ene.prss.eff.rel  = 100. * cpatch$ene.prss.eff.cum  / ene.scale
      cpatch$ene.eddy.flux.rel = 100. * cpatch$ene.eddy.flux.cum / ene.scale
      cpatch$ene.runoff.rel    = 100. * cpatch$ene.runoff.cum    / ene.scale
      cpatch$ene.residual.rel  = 100. * cpatch$ene.residual.cum  / ene.scale
      #----- Cumulative terms for all water budget. ---------------------------------------#
      cpatch$h2o.dstorage.rel  = 100. * cpatch$h2o.dstorage.cum  / h2o.scale
      cpatch$h2o.precip.rel    = 100. * cpatch$h2o.precip.cum    / h2o.scale
      cpatch$h2o.eddy.flux.rel = 100. * cpatch$h2o.eddy.flux.cum / h2o.scale
      cpatch$h2o.runoff.rel    = 100. * cpatch$h2o.runoff.cum    / h2o.scale
      cpatch$h2o.residual.rel  = 100. * cpatch$h2o.residual.cum  / h2o.scale
      #------------------------------------------------------------------------------------#







      #------------------------------------------------------------------------------------#
      #      Define a nice grid for time.                                                  #
      #------------------------------------------------------------------------------------#
      whenout = pretty.time(when,n=5)
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
         plog         = budget.now$plog
         prefix       = budget.now$prefix
         theme        = budget.now$theme
         unit         = budget.now$unit
         ylimit       = budget.now$ylim
         legpos       = budget.now$legpos
         aggrit       = budget.now$aggr
         plotit       = budget.now$plt

         if (plotit){
            #----- Define the number of layers. -------------------------------------------#
            nlayers   = length(vnames)
            namerange = vnames[lrange]
            if (any(is.na(ylimit))){
               ylimit = max(abs(cpatch[,namerange]),na.rm=TRUE)
               ylimit = c(-ylimit,ylimit)
               ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
            }#end if (any(is.na(ylimit)))
            #------------------------------------------------------------------------------#


            #----- Make plot annotations. -------------------------------------------------#
            letitre = paste0(theme," - ",thispoi$lieu
                            ,"\n Patch: ",ipa,";  Time series: ",theme)
            lex     = desc.unit(desc="Time",unit=untab$gmt    )
            ley     = desc.unit(desc=NULL  ,unit=untab[[unit]])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            cat0("        > ",theme," time series .")

            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(nout)){
               #----- Open file. ----------------------------------------------------------#
               fichier = file.path( patchdir
                                  , paste0(prefix,"-patch-",cipa,"-",suffix,".",outform[o])
                                  )#end file.path
               if(outform[o] %in% "x11"){
                  X11(width=ex.size$width,height=ex.size$height,pointsize=ptsz)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=ex.size$width*depth
                     ,height=ex.size$height*depth,pointsize=ptsz,res=depth
                     ,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=ex.size$width*depth
                      ,height=ex.size$height*depth,pointsize=ptsz,res=depth
                      ,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=ex.size$width,height=ex.size$height
                            ,pointsize=ptsz,paper=ex.size$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=ex.size$width,height=ex.size$height
                     ,pointsize=ptsz,paper=ex.size$paper)
               }#end if
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
                     , lty    = llty
                     , lwd    = llwd
                     , ncol   = 2
                     , cex    = 0.8 * cex.ptsz
                     , xpd    = TRUE
                     , bty    = "n"
                     )#end legend
               #---------------------------------------------------------------------------#


               #---- Open the window and plot the axes and annotations. -------------------#
               par(mar=c(5.1,5.1,3.1,0.6))
               plot.new()
               plot.window(xlim=range(when),ylim=ylimit,log=plog)
               axis.rt(side=1,at=whenout$levels,labels=whenout$labels,padj=whenout$padj
                      ,las=5,off=0.05)
               axis(side=2,las=1)
               if (plotgrid){
                   abline(h=axTicks(side=2),v=whenout$levels,col=grid.colour,lty="longdash")
               }#end if
               box()
               title(main=letitre,cex.main=0.6)
               title(ylab=ley,line=3.5)
               #---------------------------------------------------------------------------#



               #----- Plot the curves. ----------------------------------------------------#
               datenow = chron(ceiling(as.numeric(when)))
               today   = chron(unique(datenow))
               for (l in sequence(nlayers)){
                  #----- Aggregate data into daily groups, to reduce file size. -----------#
                  if (aggrit){
                     ynow = tapply(X=cpatch[[vnames[l]]],INDEX=datenow,FUN=mean,na.rm=TRUE)
                  }else{
                     ynow = tapply(X=cpatch[[vnames[l]]],INDEX=datenow,FUN=max ,na.rm=TRUE)
                  }#end if (aggrit)
                  #------------------------------------------------------------------------#

                  lines(x=today,y=ynow,col=lcolours[l],lwd=llwd[l],lty=llty[l]
                        ,type=ltype,pch=16)
               }#end for
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               if (outform[o] %in% "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               clean.tmp()
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
