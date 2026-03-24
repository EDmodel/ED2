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
outroot        = "thisoutroot"
#------------------------------------------------------------------------------------------#


#----- Time options. ----------------------------------------------------------------------#
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
use.trim       = myefttrim         # Use only a sub-set of years (TRUE/FALSE)
year.trima     = myeftyeara        # First year of the sub-set
year.trimz     = myeftyearz        # Last year of the sub-set
season.mona    = thisseasonmona    # First month for year by year comparison
reload.data    = TRUE              # Should I reload partially loaded data?
sasmonth.short = c(2,5,8,11)       # Months for SAS plots (short runs)
sasmonth.long  = 5                 # Months for SAS plots (long runs)
nyears.long    = 25                # Runs longer than this are considered long runs.
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
ptsz           = 16                     # Font size.
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
ibackground    = mybackground           # Background settings (check load_everything.r)
#------------------------------------------------------------------------------------------#


#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min             = -5.0         # The deepest depth that trees access water.
idbh.type           = myidbhtype   # Type of DBH class
                                   # 1 -- Every 10 cm until 100cm; > 100cm
                                   # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
klight              = myklight     # Weighting factor for maximum carbon balance
corr.growth.storage = mycorrection # Correction factor to be applied to growth and
                                   #   storage respiration
iallom              = myallom      # Allometry to use
isoil.hydro         = myslhydro    # Soil hydrology method
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


#----- Load observations. -----------------------------------------------------------------#
obsrfile = file.path(srcdir,"LBA_MIP.v9.RData")
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
   outpref = paste(outmain,"ycomp",sep="/")
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



   #----- Print a banner to entretain the user. -------------------------------------------#
   cat(" + Post-processing output from ",lieu,"...","\n")
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



   #---------------------------------------------------------------------------------------#
   #    Copy the data to temporary variables.                                              #
   #---------------------------------------------------------------------------------------#
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
   #    Create a list with unique years.                                                   #
   #---------------------------------------------------------------------------------------#
   print (paste("    - Finding the years and seasons...",sep=""))
   if (use.trim){
      sel.trim = datum$toyear %in% seq(from=year.trima,to=year.trimz,by=1)
   }else{
      sel.trim = datum$toyear > -Inf
   }#end if
   year.use     = datum$toyear[sel.trim]

   if (season.mona == 1){
      nyears    = length(year.use) + 1
      year4     = c(year.use,2*year.use[nyears-1] - year.use[nyears-2])
   }else{
      nyears    = length(year.use)
      year4     = year.use
   }#end if
   year.desc = paste(year4-c(diff(year4)[1],diff(year4)),year4,sep="-")
   year.col  = eft.col[match(year4,eft.year)]
   #----- Year for seasonal means. --------------------------------------------------------#
   yr3mon      = year.use
   nyr3mon     = length(yr3mon)
   yr3mon.desc = paste("Dec",sprintf("%2.2i"
                                    ,(yr3mon-c(diff(yr3mon)[1],diff(yr3mon))) %% 100 )
                      ,"-Nov",sprintf("%2.2i",yr3mon %% 100),sep="")
   yr3mon.col  = eft.col[match(yr3mon,eft.year)]
   yr3mon.pch  = eft.pch[match(yr3mon,eft.year)]
   #---------------------------------------------------------------------------------------#


   dwhena     = as.numeric(chron(paste(season.mona,1,year4[2],sep="/")))
   dwhenz     = as.numeric(chron(paste(season.mona,1,year4[3],sep="/")))
   sel        = datum$tomonth >= dwhena & datum$tomonth <= dwhenz
   month.when = datum$tomonth[sel]
   
   for (tt in 1:nyc.tvar){
      thisvar = yc.tvar[[tt]]
      vname   = thisvar$vnam
      desc    = thisvar$desc
      if (vname %in% c("rain","runoff","intercepted","wshed")){
         unit    = untab$mmomo
      }else{
         unit    = thisvar$unit
      }#end if
      add0             = thisvar$add
      mult0            = thisvar$mult
      colmean          = thisvar$colmean
      colerr           = thisvar$colerr
      coledge          = thisvar$coledge
      plotit           = thisvar$plt
      cumul            = thisvar$cumul
      if (cumul){
         plttype = "l"
      }else{
         plttype = "o"
      }#end if




      if (plotit){
         print(paste("    * ",desc,"..."))


         #----- Create the directories. ---------------------------------------------------#
         outtser  = paste(outpref,"tseries",sep="/")
         outbplot = paste(outpref,"bplot"  ,sep="/")
         if (! file.exists(outtser )) dir.create(outtser )
         if (! file.exists(outbplot)) dir.create(outbplot)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Create the list of data.                                                    #
         #---------------------------------------------------------------------------------#
         outplot = list()
         ylimit  = NULL
         for (y in 2:nyears){
            whena.now      = chron(paste(season.mona,1,year4[y-1],sep="/"))
            whenz.now      = chron(paste(season.mona,1,year4[y]  ,sep="/"))
            sel            = datum$tomonth >= whena.now & datum$tomonth <= whenz.now
            nsel           = sum(sel)
            outplot[[y]]   = list()
            outplot[[y]]$x = month.when
            if (cumul){
               copy.datum                    = emean[[vname]][sel]
               copy.datum[is.na(copy.datum)] = 0
               outplot[[y]]$y = cumsum(copy.datum)
            }else{
               outplot[[y]]$y = emean[[vname]][sel]
            }#end if

            #----- Update range. ----------------------------------------------------------#
            if (any(is.finite(outplot[[y]]$y)) ){
               y.min  = min(outplot[[y]]$y,na.rm=TRUE)
               y.max  = max(outplot[[y]]$y,na.rm=TRUE)
               ylimit = range(c(ylimit,y.min,y.max))
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #----- Make a dummy limit in case everything is empty. ---------------------------#
         if (is.null(ylimit)) ylimit = c(-1,1)
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #      Make the time axis.                                                        #
         #---------------------------------------------------------------------------------#
         whenplot        = pretty.time(chron(c(dwhena,dwhenz)),n=8)
         whenplot$labels = substring(mlist[nummonths(whenplot$levels)],1,3)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Plot the data.                                                             #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outtser,"/ycomp-",vname,"-",suffix,".",outform[o],sep="")
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
            
            letitre = paste(lieu," \n",desc,sep="")
            lex     = "Months"
            ley     = desc.unit(desc=desc,unit=unit)


            par(par.user)
            par(oma=c(0,0,0,0))
            layout (mat=rbind(2,1),heights=c(5,1))
            
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend(x="bottom",legend=year.desc[-1]
                  ,lwd=2.5,col=year.col[-1],title="Period",cex=0.9,ncol=4)

            par(mar=c(4.1,4.6,4.1,2.1))
            plot.new()
            plot.window(xlim=range(outplot[[2]]$x),ylim=ylimit)
            title(main=letitre,xlab=lex,ylab=ley,cex.main=0.8)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            axis(side=2,las=1)
            box()
            if (plotgrid){ 
               abline(v=whenplot$levels,h=axTicks(side=2),col=grid.colour,lty="solid")
            }#end if
            for (y in 2:nyears){
               lines(x=outplot[[y]]$x,y=outplot[[y]]$y,type=plttype,pch=16,cex=1.0
                    ,lwd=2.5,col=year.col[y])
            }#end for


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Split the data into seasonal means.                                          #
         #---------------------------------------------------------------------------------#
         print(paste("    * Seasonal bar plot: ",desc,"..."))
         this.season = season(datum$tomonth,add.year=TRUE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Make the last December become part of the first year.                        #
         #---------------------------------------------------------------------------------#
         sel              = this.season == paste(max(yr3mon)+1,"01",sep="")
         this.season[sel] = paste(min(yr3mon),"01",sep="")
         yr.season        = as.numeric(substring(this.season,1,4))
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Find the seasonality matrix.                                                #
         #---------------------------------------------------------------------------------#
         this.var = emean[[vname]]
         if (cumul){
            season.vec = tapply(X=this.var,INDEX=this.season,FUN=sum)
         }else{
            season.vec = tapply(X=this.var,INDEX=this.season,FUN=mean)
         }#end if
         season.mat = matrix( data     = season.vec
                            , ncol     = 4
                            , nrow     = nyr3mon
                            , dimnames = list(yr3mon.desc,season.list[-nseasons])
                            , byrow    = TRUE
                            )#end matrix
         if (vname %in% c("rain","runoff","intercepted","wshed")){
            ylimit = c(0,max(season.vec,na.rm=TRUE))
         }else{
            ylimit = pretty.xylim(u=season.vec,fracexp=0.0,is.log=FALSE)
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the bar plot.                                                          #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outbplot,"/bp_season_",vname,"-",suffix,".",outform[o],sep="")
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

            letitre = paste(lieu," \n Year comparison: ",desc,sep="")
            lex     = "Season"
            ley     = desc.unit(desc=desc,unit=unit)


            par(par.user)
            layout (mat=rbind(2,1),heights=c(5,1))
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend(x="bottom",inset=0.0,legend=yr3mon.desc,fill=yr3mon.col
                  ,title="Period",cex=0.9,ncol=2)


            par(mar=c(3.1,4.6,4.1,2.1))
            barplot(season.mat,col=yr3mon.col,main=letitre,xlab=lex,ylab=ley
                   ,cex.main=cex.main,ylim=ylimit,legend.text=FALSE,beside=TRUE
                   ,border=grey.fg,xpd=FALSE,las=1)
            box()


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Plot the monthly mean variables as functions of other 2 environment variables.    #
   #---------------------------------------------------------------------------------------#
   print(paste(" + Plotting parameter space: ",lieu,"...",sep=""))
   for (z in 1:nyc.zvar){
      zvname   = yc.xyzvar$zvar[[z]]$vname
      zdesc    = yc.xyzvar$zvar[[z]]$desc
      zkey     = yc.xyzvar$zvar[[z]]$key
      zunit    = yc.xyzvar$zvar[[z]]$unit

      #----- Create the directories. ------------------------------------------------------#
      outxyzp = paste(outpref,"xyzplot",sep="/")
      outzvar = paste(outxyzp,zvname   ,sep="/")
      if (! file.exists(outxyzp )) dir.create(outxyzp )
      if (! file.exists(outzvar )) dir.create(outzvar )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop over all x and y axes.                                                    #
      #------------------------------------------------------------------------------------#
      print(paste("    * Z: ",zdesc,"..."))
      for (y in 1:nyc.yvar){
         yvname   = yc.xyzvar$yvar[[y]]$vname
         ydesc    = yc.xyzvar$yvar[[y]]$desc
         yunit    = yc.xyzvar$yvar[[y]]$unit
         yleg     = yc.xyzvar$yvar[[y]]$leg

         for (x in 1:nyc.xvar){
            xvname   = yc.xyzvar$xvar[[x]]$vname
            xdesc    = yc.xyzvar$xvar[[x]]$desc
            xunit    = yc.xyzvar$xvar[[x]]$unit
            xleg     = yc.xyzvar$xvar[[x]]$leg


            print(paste("       ~ X: ",xdesc," Y: ",ydesc,"..."))

            #----- Title. -----------------------------------------------------------------#
            letitre    = paste(lieu,zdesc,sep="\n")
            #----- Attribute symbols according to the year. -------------------------------#
            this.pch  = eft.pch[match(yr.season,eft.year)]
            #----- Expand the edges of the x axis. ----------------------------------------#
            xvar      = emean[[xvname]]
            lex       = desc.unit(desc=xdesc,unit=xunit)
            xrange    = range(xvar,na.rm=TRUE)
            xlimit    = xrange
            xlimit[1] = xrange[1] - 0.05 * (xrange[2] - xrange[1])
            xlimit[2] = xrange[2] + 0.05 * (xrange[2] - xrange[1])
            #----- Expand the edges of the y axis. ----------------------------------------#
            yvar      = emean[[yvname]]
            ley       = desc.unit(desc=ydesc,unit=yunit)
            yrange    = range(yvar,na.rm=TRUE)
            ylimit    = yrange
            ylimit[1] = yrange[1] - 0.05 * (yrange[2] - yrange[1])
            ylimit[2] = yrange[2] + 0.40 * (yrange[2] - yrange[1])
            #----- Annotation for the colour map ("Z" axis). ------------------------------#
            zvar      = emean[[zvname]]
            lez       = desc.unit(desc=zdesc,unit=zunit)
            #----- Find the position to plot the legend. ----------------------------------#
            leg.pos   = paste(yleg,xleg,sep="")



            #------------------------------------------------------------------------------#
            #    Make lists for colourmap.                                                 #
            #------------------------------------------------------------------------------#
            ptitle = list()
            ptitle[[1]] = list(main=letitre,xlab=lex,ylab=ley,cex.main=cex.main)
            paxes       = list()
            paxes[[1]]  = list( x.axis = list(side=1)
                              , y.axis = list(side=2,las=1)
                              , grid   = list(col=grey.fg,lty="solid")
                              , legend = list( x      = leg.pos
                                             , inset  = 0.01
                                             , legend = yr3mon.desc
                                             , col    = foreground
                                             , bg     = background
                                             , pch    = yr3mon.pch
                                             , title  = "Period"
                                             , ncol   = 2
                                             , pt.cex = 1./0.9
                                             , cex    = 0.9
                                             )#end legend
                              )#end list
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Plot the bar plot.                                                       #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outzvar,"/cmap_x_",xvname,"_y_",yvname
                                       ,"_z_",zvname,"-",suffix,".",outform[o],sep="")
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



               #----- Plot the parameter space. -------------------------------------------#
               par(par.user)
               colourmap(x=xvar,y=yvar,z=zvar,xlim=xlimit,ylim=ylimit
                        ,colour.palette=muitas,cex=1.6,pch=this.pch,lwd=2
                        ,plot.title=ptitle
                        ,key.title=title(main=lez,cex.main=0.8)
                        ,plot.axes=paxes
                        )#end colourmap
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
               dummy = clean.tmp()
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#
