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
budget   = list()
budget[[ 1]] = list( vnam   = c("co2.dstorage","co2.nep","co2.dens.eff"
                               ,"co2.loss2atm","co2.residual")
                   , desc   = c("Delta (Storage)","NEP","Density Effect"
                               ,"Eddy flux loss","Residual")
                   , colour = c("forestgreen","chartreuse","purple4"
                               ,"deepskyblue","black")
                   , lwd    = c(2.0,2.0,2.0,2.0,2.0)
                   , range  = c(FALSE,TRUE,TRUE,TRUE,TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "carbflux"
                   , theme  = "Carbon dioxide budget"
                   , unit   = "umol/m2/s"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 2]] = list( vnam   = c("ene.dstorage","ene.precip","ene.netrad"
                               ,"ene.dens.eff","ene.prss.eff","ene.loss2atm"
                               ,"ene.drainage","ene.runoff","ene.residual")
                   , desc   = c("Delta (Storage)","Rainfall","Net Radiation"
                               ,"Density effect","Pressure effect","Eddy flux loss"
                               ,"Drainage","Runoff","Residual")
                   , colour = c("red3","royalblue","darkorange"
                               ,"purple4","chartreuse","deepskyblue","sienna"
                               ,"forestgreen","black")
                   , lwd    = c(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0)
                   , range  = c(FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "eneflux"
                   , theme  = "Enthalpy budget"
                   , unit   = "W/m2"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 3]] = list( vnam   = c("h2o.dstorage","h2o.precip","h2o.dens.eff"
                               ,"h2o.loss2atm","h2o.drainage","h2o.runoff","h2o.residual")
                   , desc   = c("Delta (Storage)","Rainfall","Density effect"
                               ,"Eddy flux loss","Drainage","Runoff","Residual")
                   , colour = c("red3","royalblue","purple4"
                               ,"deepskyblue","sienna","forestgreen","black")
                   , lwd    = c(2.0,2.0,2.0,2.0,2.0,2.0,2.0)
                   , range  = c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "h2oflux"
                   , theme  = "Water budget"
                   , unit   = "kg/m2/day"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 4]] = list( vnam   = c("co2.cumres")
                   , desc   = c("Residual")
                   , colour = c("limegreen")
                   , lwd    = c(2.0)
                   , range  = c(TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "cumco2"
                   , theme  = "CO2: cumulative residual (absolute)"
                   , unit   = "umol/m2"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 5]] = list( vnam   = c("ene.cumres")
                   , desc   = c("Residual")
                   , colour = c("red3")
                   , lwd    = c(2.0)
                   , range  = c(TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "cumene"
                   , theme  = "Enthalpy: cumulative residual (absolute)"
                   , unit   = "J/m2"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 6]] = list( vnam   = c("h2o.cumres")
                   , desc   = c("Residual")
                   , colour = c("steelblue")
                   , lwd    = c(2.0)
                   , range  = c(TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "cumh2o"
                   , theme  = "Water: cumulative residual (absolute)"
                   , unit   = "kg/m2"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 7]] = list( vnam   = c("co2.relres")
                   , desc   = c("Residual")
                   , colour = c("limegreen")
                   , lwd    = c(2.0)
                   , range  = c(TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "relco2"
                   , theme  = "CO2: cumulative residual (relative)"
                   , unit   = "---"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 8]] = list( vnam   = c("ene.relres")
                   , desc   = c("Residual")
                   , colour = c("red3")
                   , lwd    = c(2.0)
                   , range  = c(TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "relene"
                   , theme  = "Enthalpy: cumulative residual (relative)"
                   , unit   = "---"
                   , legpos = "topleft"
                   , plt    = TRUE)
budget[[ 9]] = list( vnam   = c("h2o.relres")
                   , desc   = c("Residual")
                   , colour = c("steelblue")
                   , lwd    = c(2.0)
                   , range  = c(TRUE)
                   , type   = ptype
                   , plog   = ""
                   , prefix = "relh2o"
                   , theme  = "Water: cumulative residual (relative)"
                   , unit   = "---"
                   , legpos = "topleft"
                   , plt    = TRUE)
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



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
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
   inpref     = paste(here,place,sep="/")
   outpref    = thispoi$pathout
   lieu       = thispoi$lieu
   suffix     = thispoi$iata
   #---------------------------------------------------------------------------------------#


   #----- Print the banner to entretain the user. -----------------------------------------#
   print (paste("  + ",thispoi$lieu,"...",sep=""))
   #---------------------------------------------------------------------------------------#

   #----- Make the main output directory in case it doesn't exist. ------------------------#
   if (! file.exists(outroot)) dir.create(outroot)
   outmain = paste(outroot,place,sep="/")
   if (! file.exists(outmain)) dir.create(outmain)
   outdir = paste(outmain,"budget",sep="/")
   if (! file.exists(outdir)) dir.create(outdir)
   #---------------------------------------------------------------------------------------#


   #----- Determine the number of patches. ------------------------------------------------#
   filelist  = dir(inpref)
   mypatches = length(grep("budget_state_patch_",filelist))
   #---------------------------------------------------------------------------------------#



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
      inputfile = paste(inpref,paste("budget_state_patch_",cipa,".txt",sep=""),sep="/")
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

      #----- Reduce the size of the file to be the period of interest only. ---------------#
      print(paste("      * Reduce data to the period of interest..."))
      when   = chron( chron(dates=paste(cpatch$month,cpatch$day,cpatch$year,sep="/"))
                    + cpatch$time/day.sec, out.format=c(dates="m/d/y",times="h:m:s"))
      sel    = when >= whena & when <= whenz
      cpatch = cpatch[sel,]
      when   = when[sel]

      #----- Re-scale or re-define some variables. ----------------------------------------#
      print(paste("      * Define the cumulative sum of residuals..."))
      cpatch$co2.cumres=cumsum(cpatch$co2.residual)
      cpatch$ene.cumres=cumsum(cpatch$ene.residual)
      cpatch$h2o.cumres=cumsum(cpatch$h2o.residual / day.sec)
      cpatch$co2.relres=cumsum(cpatch$co2.residual) / cpatch$co2.storage
      cpatch$ene.relres=cumsum(cpatch$ene.residual) / cpatch$ene.storage
      cpatch$h2o.relres=cumsum(cpatch$h2o.residual / day.sec) / cpatch$h2o.storage

      #------------------------------------------------------------------------------------#
      #      Define a nice grid for time.                                                  #
      #------------------------------------------------------------------------------------#
      whenout = pretty.time(when,n=8)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #   Plot the time series diagrams showing months and years.                          #
      #------------------------------------------------------------------------------------#
      print(paste("      * Plot some patch-level figures..."))
      for (bb in 1:nbudget){

         #----- Retrieve variable information from the list. ------------------------------#
         budget.now   = budget[[bb]]
         vnames       = budget.now$vnam  
         description  = budget.now$desc  
         lcolours     = budget.now$colour
         llwd         = budget.now$lwd
         lrange       = budget.now$range
         ltype        = budget.now$type
         plog         = budget.now$plog
         prefix       = budget.now$prefix
         theme        = budget.now$theme 
         unit         = budget.now$unit  
         legpos       = budget.now$legpos
         plotit       = budget.now$plt   
    
         if (plotit){


            #----- Define the number of layers. -------------------------------------------#
            nlayers   = length(vnames)
            namerange = vnames[lrange]
            ylimit = max(abs(cpatch[,namerange]),na.rm=TRUE)
            ylimit = c(-ylimit,ylimit)
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
               legend(x=legpos,inset=0.01,legend=description,col=lcolours,lwd=llwd
                     ,ncol=2,cex=0.9)
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
