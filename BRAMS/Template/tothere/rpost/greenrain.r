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
ptsz     =   15            # Font size
offnet   = 0.04            # Offset for network precipitation map.
offdcyc  = 0.25            # Offset for the diurnal cycle plots.
srcdir   = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc"
ctl2001  = "edish-leaf-_g1.ctl"  # WETAMC file
ctl2040  = "ggbare-_g1.ctl"  # WETAMC file
yearsce  = c("bogus","fixed")
seasout  = "barcab"
season   = "BARCA-B"
rainmin  = 10.
#------------------------------------------------------------------------------------------#
#    Define the sub-domain you want to use.  This is probably not going to be the exact    #
# domain, but it will be the closest it can be.                                            #
#------------------------------------------------------------------------------------------#
lona    = -79.2
lonz    = -34.8
lata    = -14.6
latz    =   9.8

plotgrid2    = FALSE
grid2.corner = list(lona=-68.0,lonz=-59.0,lata=-11.0,latz=0.0)
g2colour     = "purple4"

#----- Time to obtain the LAI, topography, and land fraction. -----------------------------#
yeara   =  2009        # Year
montha  =     5        # Month
daya    =    01        # Day
houra   =     0        # Hour

#----- Time to obtain the LAI, topography, and land fraction. -----------------------------#
yearz   =  2009        # Year
monthz  =     5        # Month
dayz    =    31        # Day
hourz   =     0        # Hour

#----- Minimum land fraction to be considered land. ---------------------------------------#
landmin =  0.9

addgrid        = FALSE     # Should I include a grid?
southamerica   = TRUE      # Is this a South America run (for map plotting)?
plotmap        = TRUE      # Should I plot the maps?
plotroad       = TRUE      # Should I plot the roads?
plotagb        = FALSE      # Should I plot the rivers?
plotriver      = FALSE     # Should I plot the rivers?
ncols          = 200       # Approximate number of colours for image plot.
colscheme      = "imuitas" # Colour scheme for precipitation rate. 
logpal         = FALSE     # Should the colour palette be in log scale? 
fontsub        = 2         # Type of font for subtitle:
                           # 1. Regular; 2. Bold; 3. Italic; 4. Bold-Italic.
cexaxes        = 1.0       # Expansion factor for axes.
keylog         = TRUE      # Should the scale be in log?
myrange        = c(10,2000)# Range to be used in the scale (NULL = default)
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

#----- Getting rid of the annoying bell. --------------------------------------------------#
options(locatorBell=FALSE)

#----- Make plotform case insensitive. ----------------------------------------------------#
plotform = tolower(plotform)

#----- Determine the number of plot output formats. ---------------------------------------#
noutput = length(plotform)
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
info01 = readctl(ctl2001)
info40 = readctl(ctl2040)

#----- Find the domain boundaries. --------------------------------------------------------#
print (" + Finding the domain boundaries...")
xa = which.min(abs(info40$glon-lona))
xz = which.min(abs(info40$glon-lonz))
ya = which.min(abs(info40$glat-lata))
yz = which.min(abs(info40$glat-latz))
ta = which.min(abs(info40$gtime-whena))
tz = which.min(abs(info40$gtime-whenz))

#----- Find the time for plot. ------------------------------------------------------------#
whena = info40$gtime[ta]
whenz = info40$gtime[tz]

#----- Read the precipitation rate. -------------------------------------------------------#
print (" + Reading BRAMS grid...")
raina01 = readgrads(vari=c("land","lai.ps","vegz0.ps","agb.ps","acccon","totpcp"),
                    info=info01,xlim=c(xa,xz),ylim=c(ya,yz),tlim=c(ta,ta))
rainz01 = readgrads(vari=c("land","lai.ps","vegz0.ps","agb.ps","acccon","totpcp"),
                    info=info01,xlim=c(xa,xz),ylim=c(ya,yz),tlim=c(tz,tz))
raina40 = readgrads(vari=c("land","lai.ps","vegz0.ps","agb.ps","acccon","totpcp"),
                    info=info40,xlim=c(xa,xz),ylim=c(ya,yz),tlim=c(ta,ta))
rainz40 = readgrads(vari=c("land","lai.ps","vegz0.ps","agb.ps","acccon","totpcp"),
                    info=info40,xlim=c(xa,xz),ylim=c(ya,yz),tlim=c(tz,tz))

limlon  = range(raina01$glon)
limlat  = range(raina01$glat)

#------------------------------------------------------------------------------------------#
#    Scale vegetation roughness to represent the rivers.                                   #
#------------------------------------------------------------------------------------------#
raina01$vegz0.ps = raina01$vegz0.ps * raina01$land + 0.0001 * (1. - raina01$land)
raina40$vegz0.ps = raina40$vegz0.ps * raina40$land + 0.0001 * (1. - raina40$land)
rainz01$vegz0.ps = rainz01$vegz0.ps * rainz01$land + 0.0001 * (1. - rainz01$land)
rainz40$vegz0.ps = rainz40$vegz0.ps * rainz40$land + 0.0001 * (1. - rainz40$land)

#------------------------------------------------------------------------------------------#
#     Organise the precipitation so it can use similar commands for absolute and relative  #
# plots.                                                                                   #
#------------------------------------------------------------------------------------------#
rainfall      = list()
rainfall[[1]] = list(lon       = raina01$glon, 
                     lat       = raina01$glat,
                     rain      = rainz01$acccon[1,1,,] + rainz01$totpcp[1,1,,]
                               - raina01$acccon[1,1,,] - raina01$totpcp[1,1,,],
                     aminusb   = FALSE,
                     maintitle = "Total precipitation",
                     period    = as.character(yearsce[1]),
                     identity  = tolower(as.character(yearsce[1])),
                     vegz0     = raina01$vegz0.ps[1,1,,],
                     land      = raina01$land[1,1,,],
                     show      = plotagb   && FALSE,
                     river     = plotriver && TRUE,
                     threshold = 1.6)
rainfall[[2]] = list(lon       = raina40$glon, 
                     lat       = raina40$glat,
                     rain      = rainz40$acccon[1,1,,] + rainz40$totpcp[1,1,,]
                               - raina40$acccon[1,1,,] - raina40$totpcp[1,1,,],
                     aminusb   = FALSE,
                     maintitle = "Total precipitation",
                     period    = as.character(yearsce[2]),
                     identity  = tolower(as.character(yearsce[2])),
                     vegz0     = raina40$vegz0.ps[1,1,,],
                     land      = raina40$land[1,1,,],
                     show      = plotagb   && TRUE,
                     river     = plotriver && FALSE,
                     threshold = 1.6)
rainfall[[3]] = list(lon       = raina40$glon, 
                     lat       = raina40$glat,
                     rain      = rainfall[[2]]$rain - rainfall[[1]]$rain,
                     maintitle = "Precipitation difference",
                     period    = tolower(paste(yearsce[2],yearsce[1],sep=" - ")),
                     identity  = "diff",
                     aminusb   = TRUE,
                     vegz0     = raina40$vegz0.ps[1,1,,],
                     land      = raina40$land[1,1,,],
                     show      = plotagb   && TRUE,
                     river     = plotriver && FALSE,
                     threshold = 1.6)
rainfall[[1]]$rain[rainfall[[1]]$rain < rainmin] = NA
rainfall[[2]]$rain[rainfall[[2]]$rain < rainmin] = NA
rainfall[[1]]$span = range(c(rainfall[[1]]$rain,rainfall[[2]]$rain),na.rm=TRUE)
rainfall[[2]]$span = rainfall[[1]]$span
rainfall[[3]]$span = c(-max(abs(rainfall[[3]]$rain),na.rm=TRUE),
                        max(abs(rainfall[[3]]$rain),na.rm=TRUE))
#-------------------------------------------------------------------------------------------#

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
      sel         = drleg > fracwhite * thisrain$span[1] & 
                    drleg < fracwhite * thisrain$span[2]
      drcleg[sel] = "white"
      keylognow   = FALSE
   }else{

      if (is.null(myrange)){
         thisspan = thisrain$span
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
      diffrainfile = paste(figpath,"/rain-",seasout,"-",thisrain$identity,".",plotform[o],
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

      letitre=paste(thisrain$maintitle," (",thisrain$period,") - ",season,sep="")
      wlaba = round(whena)
      wlabz = round(whenz)
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


                               #----- Add a grid if the user wants one. -------------------#
                               if (addgrid) grid(lty="dotted",col="lightgray")
                               #-----------------------------------------------------------#

                               #----- Add a map if the user wants one. --------------------#
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
                               #-----------------------------------------------------------#


                               #----- Add grid 2 boundary if the user wants it. -----------#
                               if (plotgrid2){
                                  rect( xleft   = grid2.corner$lona
                                      , ybottom = grid2.corner$lata
                                      , xright  = grid2.corner$lonz
                                      , ytop    = grid2.corner$latz
                                      , lwd = 2, border=g2colour)
                               }#end if

                               #-----------------------------------------------------------#
                               #     Plot the AGB to define deforestation if this is       #
                               # sought by the user.                                       #
                               #-----------------------------------------------------------#
                               if (thisrain$show){
                                  contour(x=thisrain$lon,y=thisrain$lat,z=thisrain$vegz0,
                                          levels=thisrain$threshold,drawlabels=FALSE,
                                          col=agbopt$col,lty=agbopt$lty,add=TRUE,
                                          lwd=agbopt$lwd)
                               }#end if
                               #-----------------------------------------------------------#


                               #----- Add the contours to tell where the forest is. -------#
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

                               #----- Add some location info. -----------------------------#
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
