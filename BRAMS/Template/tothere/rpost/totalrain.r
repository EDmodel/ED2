removeall = TRUE
if (removeall){
   print ("Removing...")
   rm(list=ls())
   graphics.off()
   recomp = TRUE
}else{
   graphics.off()
   recomp = FALSE
}#end if

graphics.off()

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#                                  CHANGE THE SETTINGS HERE.                               #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#----- Output options. --------------------------------------------------------------------#
path    = getwd()          # Working directory.
plotform = c("png")        # Type of output: currently valid options are png, eps, and x11.
paper    = "letter"        # Paper size
depth    =   96            # Resolution in pixels per inch (png only)
ptsz     =   12            # Font size
offnet   = 0.04            # Offset for network precipitation map.
offdcyc  = 0.25            # Offset for the diurnal cycle plots.
srcdir   = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc"
ctlfiles = c("abarca-_g1.ctl","bbarca-_g1.ctl")
season   = c("BARCA A","BARCA B")
slabel   = c("abarca","bbarca")
rainmin  = 10.
#------------------------------------------------------------------------------------------#
#    Define the sub-domain you want to use.  This is probably not going to be the exact    #
# domain, but it will be the closest it can be.                                            #
#------------------------------------------------------------------------------------------#
lona    = -79.2
lonz    = -34.8
lata    = -14.6
latz    =   9.8

#----- Time domain, for both periods. -----------------------------------------------------#
yeara   =  c(2008, 2009)        # Year
montha  =  c(  11,    5)        # Month
daya    =  c(  14,    7)        # Day
houra   =  c(   0,    0)        # Hour

#----- Time to obtain the LAI, topography, and land fraction. -----------------------------#
yearz   =  c(2008, 2009)        # Year
monthz  =  c(  12,    5)        # Month
dayz    =  c(   7,   31)        # Day
hourz   =  c(   0,    0)        # Hour


#----- Minimum land fraction to be considered land. ---------------------------------------#
landmin =  0.9

addgrid        = FALSE     # Should I include a grid?
southamerica   = TRUE      # Is this a South America run (for map plotting)?
plotmap        = TRUE      # Should I plot the maps?
plotroad       = TRUE      # Should I plot the roads?
plotagb        = FALSE     # Should I plot the rivers?
plotriver      = FALSE     # Should I plot the rivers?
ncols          = 200       # Approximate number of colours for image plot.
colscheme      = "imuitas" # Colour scheme for precipitation rate. 
logpal         = TRUE      # Should the colour palette be in log scale? 
fontsub        = 2         # Type of font for subtitle:
                           # 1. Regular; 2. Bold; 3. Italic; 4. Bold-Italic.
cexaxes        = 1.0       # Expansion factor for axes.
keylog         = TRUE      # Should the scale be in log?
myrange        = NULL      # Range to be used in the scale (NULL = default)
fracwhite      = 0.10      # Fraction too close to zero to be plotted

#----- List of points to be plotted. -------------------------------------------------------#
site = list()
site[[1]] = list(iata="PVH",lon=-63.90,lat=-8.71)
site[[2]] = list(iata="MAO",lon=-60.04,lat=-3.03)
site[[3]] = list(iata="STM",lon=-54.71,lat=-2.44)
site[[4]] = list(iata="BEL",lon=-48.50,lat=-1.46)
site[[5]] = list(iata="AFL",lon=-56.09,lat=-9.88)
site[[6]] = list(iata="BVB",lon=-60.65,lat= 2.82)
site[[7]] = list(iata="SJL",lon=-67.09,lat=-0.13)

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#                        UNLESS YOU ARE DEVELOPING THIS SCRIPT...                          #
#                 YOU DON'T NEED TO MODIFY ANYTHING BEYOND THIS POINT!!!                   #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#----- Load libraries. --------------------------------------------------------------------#
isok = require(chron)
isok = require(fields)
isok = require(akima)

#----- Define useful paths. ---------------------------------------------------------------#
mapdir = paste(srcdir,"samap",sep="/")
rpath  = paste(srcdir,"roads",sep="/")

#----- Load functions in Rsc. -------------------------------------------------------------#
source(paste(srcdir,"atlas.r",sep="/"))
source(paste(srcdir,"bintemp.r",sep="/"))
source(paste(srcdir,"eps.r",sep="/"))
source(paste(srcdir,"gridp.r",sep="/"))
source(paste(srcdir,"gridt.r",sep="/"))
source(paste(srcdir,"muitas.r",sep="/"))
source(paste(srcdir,"peb.r",sep="/"))
source(paste(srcdir,"plotsize.r",sep="/"))
source(paste(srcdir,"pretty.log.r",sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
source(paste(srcdir,"readctl.r",sep="/"))
source(paste(srcdir,"readgrads.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"thermlib.r",sep="/"))
source(paste(srcdir,"roadmap.r",sep="/"))

#----- The list of variables to be loaded from the file. ----------------------------------#
vlist = c("land","lai.ps","vegz0.ps","agb.ps","acccon","totpcp")

#----- Getting rid of the annoying bell. --------------------------------------------------#
options(locatorBell=FALSE)

#----- Make plotform case insensitive. ----------------------------------------------------#
plotform = tolower(plotform)

#----- Determine the number of plot output formats. ---------------------------------------#
noutput = length(plotform)
nctl    = length(ctlfiles)
nsites  = length(site)

#----- Find the time for plot. ------------------------------------------------------------#
whena = chron(dates=paste(montha,daya,yeara,sep="/"),times=paste(houra,0,0,sep=":"))
whenz = chron(dates=paste(monthz,dayz,yearz,sep="/"),times=paste(hourz,0,0,sep=":"))

#------------------------------------------------------------------------------------------#
#      Determine the output and plot output directories, and check whether they exist or   #
# not.  In case they don't, create them.                                                   #                
#------------------------------------------------------------------------------------------#
print (" + Checking output paths...")
figpath = paste(path,"figures",sep="/")
if (! file.exists(figpath)) dir.create(figpath)

#----- Retrieving the control file. -------------------------------------------------------#
print (" + Reading BRAMS control files...")
info    = list()
for (n in 1:nctl) info[[n]] = readctl(ctlfiles[n])

#----- Find the domain boundaries. --------------------------------------------------------#
print (" + Finding the domain boundaries...")
xa = which.min(abs(info[[1]]$glon-lona))
xz = which.min(abs(info[[1]]$glon-lonz))
ya = which.min(abs(info[[1]]$glat-lata))
yz = which.min(abs(info[[1]]$glat-latz))
ta = tz = c(NA,NA)
for (n in 1:nctl){
  ta[n] = which.min(abs(info[[n]]$gtime-whena[n]))
  tz[n] = which.min(abs(info[[n]]$gtime-whenz[n]))
}#end for

#----- Find the time for plot. ------------------------------------------------------------#
for (n in 1:nctl){
   whena[n] = info[[n]]$gtime[ta[n]]
   whenz[n] = info[[n]]$gtime[tz[n]]
}#end if

limlon  = range(info[[1]]$glon[xa:xz])
limlat  = range(info[[1]]$glat[ya:yz])


if (recomp){
   #----- Initialise the list. ------------------------------------------------------------#
   data       = list()
   aver       = list()
   namesaver  = list()
   dmean      = list()
   namesdmean = list()
   #----- Read the data from all simulations. ---------------------------------------------#
   print (" + Reading BRAMS grid...")
   for (n in 1:nctl){
      data[[n]] = list()
      data[[n]][[1]]  = readgrads(vari=vlist,info=info[[n]],xlim=c(xa,xz),ylim=c(ya,yz),
                                  tlim=c(ta[n],ta[n]))
      data[[n]][[2]]  = readgrads(vari=vlist,info=info[[n]],xlim=c(xa,xz),ylim=c(ya,yz),
                                  tlim=c(tz[n],tz[n]))

      for (x in 1:length(data[[n]])){
         data[[n]][[x]]$vegz0.ps = data[[n]][[x]]$vegz0.ps * data[[n]][[x]]$land +
                                   0.0001                  * (1. - data[[n]][[x]]$land)
      }#end for
   }#end for

   #---------------------------------------------------------------------------------------#
   #     Organise the precipitation so it can use similar commands for absolute and rel-   #
   # ative plots.                                                                          #
   #---------------------------------------------------------------------------------------#
   rainfall      = list()
   rainspan      = NULL
   for (n in 1:nctl){
      rainfall[[n]] = list(lon       = data[[n]][[1]]$glon, 
                           lat       = data[[n]][[1]]$glat,
                           rain      = data[[n]][[2]]$acccon[1,1,,] + 
                                       data[[n]][[2]]$totpcp[1,1,,] - 
                                       data[[n]][[1]]$acccon[1,1,,]  - 
                                       data[[n]][[1]]$totpcp[1,1,,],
                           aminusb   = FALSE,
                           maintitle = "Total precipitation",
                           period    = as.character(season[n]),
                           identity  = as.character(slabel[n]),
                           vegz0     = data[[n]][[1]]$vegz0.ps[1,1,,],
                           land      = data[[n]][[1]]$land[1,1,,],
                           show      = plotagb,
                           river     = plotriver,
                           threshold = 1.6)
      rainfall[[n]]$rain[rainfall[[n]]$rain < rainmin] = NA
      rainspan = range(c(rainspan,rainfall[[n]]$rain),na.rm=TRUE)
   }#end for
}#end if
nrain=length(rainfall)

#------------------------------------------------------------------------------------------#
#      Plotting section begins here...                                                     #
#------------------------------------------------------------------------------------------#
print ("    - Defining the colours of the maps...")
#----- Define some variables based on the combination of layers. --------------------------#
if (plotmap){
   if (plotagb && plotroad){
      mapopt  = list(lty = "dotted", col = "black" , lwd = 2)
      agbopt  = list(lty = "solid" , col = "gray29", lwd = 2)
      roadopt = list(lty = "solid" , col = "black" , lwd = 2)
   }else if (plotagb){
      mapopt  = list(lty = "solid" , col = "gray29", lwd = 2)
      agbopt  = list(lty = "solid" , col = "black" , lwd = 2)
      roadopt = list(lty = "solid" , col = "white" , lwd = 2)
   }else if (plotroad){
      mapopt  = list(lty = "solid" , col = "gray29", lwd = 2)
      agbopt  = list(lty = "solid" , col = "white" , lwd = 2)
      roadopt = list(lty = "solid" , col = "black" , lwd = 2)
   }else{
      mapopt  = list(lty = "solid" , col = "black" , lwd = 2)
      agbopt  = list(lty = "solid" , col = "white" , lwd = 2)
      roadopt = list(lty = "solid" , col = "white" , lwd = 2)
   }#end if
}else if(plotagb && plotroad){
   mapopt  = list(lty = "solid" , col = "white" , lwd = 2)
   agbopt  = list(lty = "solid" , col = "gray29", lwd = 2)
   roadopt = list(lty = "solid" , col = "black" , lwd = 2)
}else if (plotagb){
   mapopt  = list(lty = "solid" , col = "white" , lwd = 2)
   agbopt  = list(lty = "solid" , col = "black" , lwd = 2)
   roadopt = list(lty = "solid" , col = "white" , lwd = 2)
}else if (plotroad){
   mapopt  = list(lty = "solid" , col = "white" , lwd = 2)
   agbopt  = list(lty = "solid" , col = "white" , lwd = 2)
   roadopt = list(lty = "solid" , col = "black" , lwd = 2)
}else{
   mapopt  = list(lty = "solid" , col = "white" , lwd = 2)
   agbopt  = list(lty = "solid" , col = "white" , lwd = 2)
   roadopt = list(lty = "solid" , col = "white" , lwd = 2)
}#end if
#------------------------------------------------------------------------------------------#


#----- Read the precipitation rate. --------------------------------------------------------#
print (" + Defining the plot size...")
size = plotsize(proje=TRUE,limlon=limlon,limlat=limlat,extendfc=TRUE,paper=paper)

for (rr in 1:nrain){
   thisrain = rainfall[[rr]]

   if (thisrain$aminusb){
      drleg       = pretty(thisrain$span,n=ncols)
      drcleg      = imuitas(n=length(drleg))
      sel         = drleg > fracwhite * drainspan[1] & 
                    drleg < fracwhite * drainspan[2]
      drcleg[sel] = "white"
      keylognow   = FALSE
   }else{

      if (is.null(myrange)){
         thisspan = rainspan
      }else{
         thisspan = myrange
      }#end if

      if (keylog){
         drleg       = pretty.log(thisspan,n=ncols)
      }else{
         drleg       = pretty(thisspan,n=ncols)
      }
      drcleg      = imuitas(n=length(drleg))
      keylognow   = keylog
   }#end if
   ndr = length(drcleg)


   for (o in 1:noutput){
      diffrainfile = paste(figpath,"/rain-",thisrain$identity,".",plotform[o],
                           sep="")
      if (plotform[o] == "x11"){
         X11(width=size$width,height=size$height,pointsize=ptsz)
      }else if(plotform[o] == "eps"){
          postscript(file=diffrainfile,width=size$width,height=size$height,
                     pointsize=ptsz,paper=paper)
      }else{
          png(filename=diffrainfile,width=size$width*depth,height=size$height*depth,
              res=depth,pointsize=ptsz)
      }#end if

      letitre=paste(thisrain$maintitle," (",thisrain$period,")",sep="")
      wlaba = round(whena[rr])
      wlabz = round(whenz[rr])
      wlaba = chron(wlaba,out.format=c(dates="mon day"))
      wlabz = chron(wlabz,out.format=c(dates="mon day"))
      lesub=paste("Period: ",wlaba," - ",wlabz,sep="")
      sombreado(x=thisrain$lon,y=thisrain$lat,z=thisrain$rain,
                levels=drleg,nlevels=ndr,col=drcleg,
                plot.title = title(main=letitre,sub=lesub,xlab="",ylab="",font.sub=2),
                key.title  = title(main="[mm]",cex.main=0.8),
                key.axes   = axis(side=4,cex.axis=cexaxes),
                key.log    = keylognow,
                plot.axes  = { axis(side=1,cex.axis=cexaxes)
                               axis(side=2,cex.axis=cexaxes)


                               #----- Add a grid if the user wants one. -----------------#
                               if (addgrid) grid(lty="dotted",col="lightgray")
                               #---------------------------------------------------------#

                               #----- Add a map if the user wants one. ------------------#
                               if (plotmap){
                                  if (southamerica){
                                     southammap(mapdir=mapdir,col=mapopt$col
                                               ,lty=mapopt$lty,lwd=mapopt$lwd)
                                  }else{
                                     map("state",fill=FALSE,boundary=TRUE,add=TRUE
                                        ,col=mapopt$col,lty=mapopt$lty,lwd=mapopt$lwd)
                                     map("worldHires",fill=FALSE,boundary=TRUE,add=TRUE
                                        ,col=mapopt$col,lty=mapopt$lty,lwd=mapopt$lwd)
                                  }#end if
                               }#end if
                               #---------------------------------------------------------#

                               #---------------------------------------------------------# 
                               #     Plot the AGB to define deforestation if this is     #
                               # sought by the user.                                     #
                               #---------------------------------------------------------#
                               if (thisrain$show){
                                  contour(x=thisrain$lon,y=thisrain$lat,z=thisrain$vegz0,
                                          levels=thisrain$threshold,drawlabels=FALSE,
                                          col=agbopt$col,lty=agbopt$lty,add=TRUE,
                                          lwd=agbopt$lwd)
                               }#end if
                               #-----------------------------------------------------------#


                               #----- Add the contours to tell where the forest is. --------#
                               if (thisrain$river){
                                  contour(x=thisrain$lon,y=thisrain$lat,z=thisrain$land,
                                          levels=landmin,drawlabels=FALSE,
                                          col="darkviolet",add=TRUE,lwd=2,lty="dotted")
                               }#end if


                               #----- Add the road map if the user wants it. --------------#
                               if (plotroad){
                                  roadmap(rpath=rpath,lty=roadopt$lty,lwd=roadopt$lwd
                                         ,col=roadopt$col)
                               }#end if
                               #-----------------------------------------------------------#

                               #----- Add some location info. ------------------------------#
                               for (ss in 1:nsites){
                                  thissite = site[[ss]]
                                  text(x=thissite$lon,y=thissite$lat,labels=thissite$iata
                                      ,cex=1.0,font=2,col="gray16",adj=c(0.5,0.5),srt=0)
                               }#end for
                          } #end plot axes
                ) #end sombreado

      #----- Close device. -------------------------------------------------------------------#
      if (plotform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
   }#end for
}#end for
