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
srcdir   = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc"
ctl      = c("abarca-_g1.ctl","bbarca-_g1.ctl")
simpref  = c("abarca","bbarca")
simdescr = c("BARCA-A","BARCA-B")
#----- List of variables to be plotted. ---------------------------------------------------#
vardf        = data.frame(t(rep(NA,times=6)))
names(vardf) = c("name","outpref","description","unit","colour","smallfrac")
vardf[ 1,] = c("rlong"   ,"rlong"   ,"Downward longwave Radiation"  ,     "W/m2"
              ,"cloudy",0.10)
vardf[ 2,] = c("qwflxca" ,"latent"  ,"Evapotranspiration"           ,  "mm/day" 
              ,"imuitas" ,0.10)
vardf[ 3,] = c("pblhgt" ,"pblhgt"   ,"PBL Height"                   ,  "     m" 
              ,"muitas"  ,0.10)
#------------------------------------------------------------------------------------------#



#----- List of points to be plotted. ------------------------------------------------------#
site = list()
site[[ 1]] = list(iata="PVH",lon=-63.90,lat=-8.71)
site[[ 2]] = list(iata="MAO",lon=-60.04,lat=-3.03)
site[[ 3]] = list(iata="STM",lon=-54.71,lat=-2.44)
site[[ 4]] = list(iata="BEL",lon=-48.50,lat=-1.46)
site[[ 5]] = list(iata="AFL",lon=-56.09,lat=-9.88)
site[[ 6]] = list(iata="BVB",lon=-60.65,lat= 2.82)
site[[ 7]] = list(iata="SJL",lon=-67.09,lat=-0.13)
site[[ 8]] = list(iata="TFF",lon=-64.72,lat=-3.38)
site[[ 9]] = list(iata="RBR",lon=-67.80,lat=-9.99)
site[[10]] = list(iata="HIA",lon=-63.07,lat=-7.53)
site[[11]] = list(iata="LBR",lon=-64.77,lat=-7.28)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Define the sub-domain you want to use.  This is probably not going to be the exact    #
# domain, but it will be the closest it can be.                                            #
#------------------------------------------------------------------------------------------#
lona    = -79.2
lonz    = -34.8
lata    = -14.6
latz    =   9.8
nzg     =    12

#----- Time to obtain the LAI, topography, and land fraction. -----------------------------#
yeara   = c( 2008, 2009)        # Year
montha  = c(   11,    5)        # Month
daya    = c(   14,    7)        # Day
houra   = c(    0,    0)        # Hour

#----- Time to obtain the LAI, topography, and land fraction. -----------------------------#
yearz   = c( 2008, 2009)        # Year
monthz  = c(   12,    5)        # Month
dayz    = c(    7,   31)        # Day
hourz   = c(    0,    0)        # Hour

#----- Minimum land fraction to be considered land. ---------------------------------------#
landmin =  0.9

addgrid        = FALSE     # Should I include a grid?
hovgrid        = TRUE      # Should I plot a grid in the POI plots?
southamerica   = TRUE      # Is this a South America run (for map plotting)?
plotmap        = TRUE      # Should I include a map?
plotriver      = ! plotmap # Should I include rivers?
plotroad       = TRUE      # Should I include the roads?
plotdef        = TRUE      # Should I include the deforestation.
ncolours       = 200       # Approximate number of colours for image plot.
colscheme      = "imuitas" # Colour scheme for precipitation rate. 
logpal         = FALSE     # Should the colour palette be in log scale? 
fontsub        = 2         # Type of font for subtitle:

defthreshold   = 1.6       # Threshold for disturbed vegetation (roughness)
cexaxes        = 1.0       # Expansion factor for axes.

myplaces       =  c("pvelho","manaus","tefe","santarem","sao_gabriel","belem","sinop"
                   ,"serra_do_navio")

#----- Similar to Hovmoller diagrams. -----------------------------------------------------#
nhov = 0
hovdi01 = list(vnam   = c("gpp","plresp","resphet","nep")
              ,desc   = c("GPP","Plant resp.","Het. resp.","NEP")
              ,colour = c("forestgreen","chartreuse","sienna","darkolivegreen")
              ,prefix = "carbflux"
              ,theme  = "Ecosystem carbon fluxes"
              ,unit   = "kgC/m2/yr"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi02 = list(vnam   = c("co2","xo2can.ps")
              ,desc   = c("Atmosphere","Canopy air")
              ,colour = c("deepskyblue","gray21")
              ,prefix = "co2"
              ,theme  = "Carbon dioxide mixing ratio"
              ,unit   = "ppm"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi03 = list(vnam   = c("tempc","tcan.ps","tveg.ps")
              ,desc   = c("Atmosphere","Canopy air","Leaf")
              ,colour = c("deepskyblue","gray21","chartreuse")
              ,prefix = "temperature"
              ,theme  = "Temperature"
              ,unit   = "degC"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi04 = list(vnam   = c("rv","rvcan.ps")
              ,desc   = c("Atmosphere","Canopy air")
              ,colour = c("deepskyblue","gray21")
              ,prefix = "h2ovapour"
              ,theme  = "Water vapour mixing ratio"
              ,unit   = "g/kg"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi05 = list(vnam   = c("rv","rvcan.ps")
              ,desc   = c("Atmosphere","Canopy air")
              ,colour = c("deepskyblue","gray21")
              ,prefix = "h2ovapour"
              ,theme  = "Water vapour mixing ratio"
              ,unit   = "g/kg"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi06 = list(vnam   = c("rsnet","rlong","rlongup","hflxca","qwflxca")
              ,desc   = c("Net SW radiation","Down LW radiation","Up LW radiation"
                         ,"Sensible heat","Latent heat")
              ,colour = c("darkorange","deepskyblue","sienna","midnightblue","chartreuse")
              ,prefix = "enflux"
              ,theme  = "Energy fluxes"
              ,unit   = "W/m2"
              ,legpos = "topleft"
              ,plt    = TRUE)

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

#----- List of variables. -----------------------------------------------------------------#
nvars  = length(vardf[,1])
nctl   = length(ctl)
nlist  = nctl + 1
nsites = length(site)
npois  = length(myplaces)

#----- Load functions in Rsc. -------------------------------------------------------------#
source(paste(srcdir,"atlas.r",sep="/"))
source(paste(srcdir,"bintemp.r",sep="/"))
source(paste(srcdir,"cloudy.r",sep="/"))
source(paste(srcdir,"eps.r",sep="/"))
source(paste(srcdir,"globdims.r",sep="/"))
source(paste(srcdir,"gridp.r",sep="/"))
source(paste(srcdir,"gridt.r",sep="/"))
source(paste(srcdir,"locations.r",sep="/"))
source(paste(srcdir,"muitas.r",sep="/"))
source(paste(srcdir,"peb.r",sep="/"))
source(paste(srcdir,"plotsize.r",sep="/"))
source(paste(srcdir,"pretty.log.r",sep="/"))
source(paste(srcdir,"qapply.r",sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
source(paste(srcdir,"readctl.r",sep="/"))
source(paste(srcdir,"readgrads.r",sep="/"))
source(paste(srcdir,"roadmap.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"thermlib.r",sep="/"))
source(paste(srcdir,"timeutils.r",sep="/"))

#----- Getting rid of the annoying bell. --------------------------------------------------#
options(locatorBell=FALSE)

#----- Make plotform case insensitive. ----------------------------------------------------#
plotform = tolower(plotform)

#----- Determine the number of plot output formats. ---------------------------------------#
noutput = length(plotform)

#----- Find the time for plot. ------------------------------------------------------------#
whena = chron(dates=paste(montha,daya,yeara,sep="/"),times=paste(houra,0,0,sep=":"))
whenz = chron(dates=paste(monthz,dayz,yearz,sep="/"),times=paste(hourz,0,0,sep=":"))

#----- Multiple time-series diagram -------------------------------------------------------#
if (nhov > 0){
    hovdi      = list()
    nameshovdi = NULL
    for (s in 1:nhov){
      sss        = substring(100+s,2,3)
      hdhd       = paste("hovdi",sss,sep="")
      nameshovdi = c(nameshovdi,hdhd)
      hovdi[[s]] = get(hdhd)
    } #end for
    names(hovdi) = nameshovdi
}#end if
#------------------------------------------------------------------------------------------#

#----- The list of variables to be loaded from the file. ----------------------------------#
# vlist = c("lon","lat","lai.ps","agb.ps","vegz0.ps","land","ue.avg","ve.avg"
#          ,"hflxca","qwflxca","rshort","albedt","rlong","rlongup","gpp"
#          ,"plresp","resphet","tcan.ps","rvcan.ps","co2can.ps","co2","rv"
#          ,"tempc","tsoil.ps","tveg.ps","conprr01","pcprate")
vlist = c("lon","lat","lai.ps","agb.ps","vegz0.ps","land","qwflxca","pblhgt","rlong")

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
for (n in 1:nctl) info[[n]] = readctl(ctl[n])

#----- Find the domain boundaries. --------------------------------------------------------#
print (" + Finding the domain boundaries...")
xa = which.min(abs(info[[1]]$glon-lona))
xz = which.min(abs(info[[1]]$glon-lonz))
ya = which.min(abs(info[[1]]$glat-lata))
yz = which.min(abs(info[[1]]$glat-latz))
ta = tz = rep(NA,times=nctl)
for (n in 1:nctl){
   ta[n] = which.min(abs(info[[n]]$gtime-whena[n]))
   tz[n] = which.min(abs(info[[n]]$gtime-whenz[n]))
}#end for

limlon  = range(info[[1]]$glon[xa:xz])
limlat  = range(info[[1]]$glat[ya:yz])

#----- Read the precipitation rate. -------------------------------------------------------#
print (" + Defining the plot size...")
size = plotsize(proje=TRUE,limlon=limlon,limlat=limlat,extendfc=TRUE,paper=paper)


if (recomp){
   #----- Initialise the list. ------------------------------------------------------------#
   data       = list()
   aver       = list()
   dmean      = list()

   #----- Read the data from all simulations. ---------------------------------------------#
   print (" + Reading BRAMS grid...")
   first = TRUE
   for (n in 1:nctl){
    
      dd = readgrads(vari=vlist,info=info[[n]],xlim=c(xa,xz),ylim=c(ya,yz),zlim=c(1,nzg),
                     tlim=c(ta[n],tz[n]))
      dd$zmax[dd$zmax > 1] = 1

      if (first){
         print ("   - Retrieving the nearest neighbour of the following POIs...")
         poi = list()
         for (ipy in 1:npois){
            place = myplaces[ipy]
            thispoi  = locations(where=place,here=path,monthly=TRUE,fullonly=FALSE)

            pp          = list()
            pp$pathroot = paste(path,"figures",place,sep="/")
            pp$outpref  = pp$pathroot
            pp$lieu     = thispoi$lieu
            pp$suffix   = thispoi$iata
            pp$lon      = thispoi$lon
            pp$lat      = thispoi$lat

            print (paste("     * ",pp$lieu,"...",sep=""))
            if (! file.exists(pp$outpref))  dir.create(pp$outpref)

            alllon = as.vector(dd$lon[1,1,,])
            alllat = as.vector(dd$lat[1,1,,])
            xcoord = rep(seq(from=1,to=dd$xmax),times=dd$ymax)
            ycoord = rep(seq(from=1,to=dd$ymax),each=dd$xmax)
            icl    = which.min(rdist.earth(x1=cbind(alllon,alllat)
                                          ,x2=cbind(pp$lon,pp$lat)
                                          ,miles=FALSE))
            pp$xatm     = xcoord[icl]
            pp$yatm     = ycoord[icl]

            poi[[ipy]]          = pp
         }#end for

         first = FALSE
      }#end if


      #----- Find the hours. --------------------------------------------------------------#
      dd$hour  = hours(dd$gtime)
      #------------------------------------------------------------------------------------#

      #----- Read the data from all simulations. ------------------------------------------#
      allhrs          = sort(unique(dd$hour))
      data[[n]]       = dd
      aver[[n]]       = list(xmax=dd$xmax,ymax=dd$ymax,zmax=dd$zmax
                            ,glev=dd$glev,glon=dd$glon,glat=dd$glat)
      dmean[[n]]      = list(xmax=dd$xmax,ymax=dd$ymax,zmax=dd$zmax
                            ,tmax=length(allhrs),hour=allhrs,glev=dd$glev
                            ,glon=dd$glon,glat=dd$glat)
      #------------------------------------------------------------------------------------#
   }#end for

   #----- Initialise the lists which will contain the differences. ------------------------#
   print (" + Computing the instantaneous differences between scenarios...")
   data[[nlist]]       = list(tmax=data[[nctl]]$tmax,zmax=data[[nctl]]$zmax
                             ,xmax=data[[nctl]]$xmax,ymax=data[[nctl]]$ymax
                             ,gtime=data[[nctl]]$gtime,glev=data[[nctl]]$glev
                             ,glon=data[[nctl]]$glon,glat=data[[nctl]]$glat
                             ,hour=data[[nctl]]$hour
                             ,lon=data[[nctl]]$lon,lat=data[[nctl]]$lat
                             ,lai.ps=data[[nctl]]$lai.ps,agb.ps=data[[nctl]]$agb.ps
                             ,vegz0.ps=data[[nctl]]$vegz0.ps,land=data[[nctl]]$land)
   aver[[nlist]]       = aver[[1]]
   dmean[[nlist]]      = dmean[[1]]
   #----- Loop over all variables to compute the diurnal cycle. ---------------------------#
   print(" + Computing the seasonal averages...")
   for (vv in 1:nvars){
      thisvar   = vardf$name[vv]
      thisdescr = vardf$description[vv]
      print (paste("   - Variable: ",thisdescr,"...",sep=""))

      #----- Copy the variable to a scratch matrix. ---------------------------------------#
      aux      = list()
      for (n in 1:nctl){
         aux[[n]] = matrix(data[[n]][[thisvar]],nrow=data[[n]]$tmax,
                                                ncol=data[[n]]$xmax*data[[n]]$ymax)
         seasmean = apply(X=aux[[n]]  ,MARGIN=2,FUN=mean,na.rm=TRUE)
         dcycmean = qapply(mat=aux[[n]],bycol=TRUE,index=data[[n]]$hour
                          ,func=mean,na.rm=TRUE)

         #----- Compute the mean differences, then bind the differences to the list. ------#
         aver[[n]][[thisvar]]  = matrix(seasmean,nrow=aver[[n]]$xmax  ,ncol=aver[[n]]$ymax)
         dmean[[n]][[thisvar]] = array(dcycmean,dim=c(dmean[[n]]$tmax,dmean[[n]]$xmax,
                                                      dmean[[n]]$ymax))
      }#end for
      aver[[nlist]][[thisvar]]  = aver[[2]][[thisvar]]  - aver[[1]][[thisvar]]
      dmean[[nlist]][[thisvar]] = dmean[[2]][[thisvar]] - dmean[[1]][[thisvar]]
   }#end for


   #----- Loop over all polygons to copy the diurnal cycle. -------------------------------#
   print(" + Copying the POIs' diurnal cycle...")
   for (ipy in 1:npois){
      pp       = poi[[ipy]]
      print(paste("   - ",pp$lieu,"...",sep=""))

      pp$tmax  = dmean[[1]]$tmax
      pp$hour  = dmean[[1]]$hour
      pp$dmean = list()
      for (n in 1:nlist){
         pp$dmean[[n]] = list()
         for (vv in 1:nvars){
            thisvar   = vardf$name[vv]
            thisdescr = vardf$description[vv]
            pp$dmean[[n]][[thisvar]] = dmean[[n]][[thisvar]][,pp$xatm,pp$yatm]
         }#end for
      }#end for
      poi[[ipy]] = pp
   }#end for
}#end if




#----- Find the range for each variable. --------------------------------------------------#
rangevar           = matrix(NA,ncol=6,nrow=nvars)
dimnames(rangevar) = list(vardf$name,c("smin","smax","smaxabs","dmin","dmax","dmaxabs"))
for (vv in 1:nvars){
   thisvar        = vardf$name[vv]
   allaver        = c(aver[[1]][[thisvar]],aver[[2]][[thisvar]])
   alldmean       = c(dmean[[1]][[thisvar]],dmean[[2]][[thisvar]])
   rangevar[vv,1] = min(allaver                   ,na.rm=TRUE)
   rangevar[vv,2] = max(allaver                   ,na.rm=TRUE)
   rangevar[vv,3] = max(abs(aver[[3]][[thisvar]]) ,na.rm=TRUE)
   rangevar[vv,4] = min(alldmean                  ,na.rm=TRUE)
   rangevar[vv,5] = max(alldmean                  ,na.rm=TRUE)
   rangevar[vv,6] = max(abs(dmean[[3]][[thisvar]]),na.rm=TRUE)
}#end for
#------------------------------------------------------------------------------------------#



#----- Define some instructions on what and how to plot the averages. ---------------------#
titlelist    = c("Seasonal mean","Seasonal mean","Mean difference")
plotlai      = c(FALSE, TRUE,TRUE)
difference   = c(FALSE,FALSE,TRUE)
plotdefnow   = plotdef & c(FALSE,TRUE,TRUE)
plotrivernow = plotriver & c(TRUE,FALSE,FALSE)
seasonwhen   = c(simdescr,paste(simdescr[2],simdescr[1],sep=" - "))
suffout      = c(simpref,"diff")
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Create the output for the three averages, and all variables.                          #
#------------------------------------------------------------------------------------------#
for (aa in 1:nlist){

   print (paste(" + Plotting full-time seasonal average for ",seasonwhen[aa],"...",sep=""))

   #----- Retrieve information about this average. ----------------------------------------#
   avenow     = aver[[aa]]
   titlenow   = titlelist[aa]
   plotlainow = plotlai[aa]
   diffnow    = difference[aa]
   suffoutnow = suffout[aa]
   seasonnow  = seasonwhen[aa]

   for (vv in 1:nvars){

      #---- Find the output-related variable information. ---------------------------------#
      thisvar     = vardf$name        [vv]
      thisoutpref = vardf$outpref     [vv]
      thisdescr   = vardf$description [vv]
      thisunit    = vardf$unit        [vv]
      thiscolour  = get(vardf$colour  [vv])
      thissmfrac  = as.numeric(vardf$smallfrac   [vv])


      print (paste("   - Variable: ",thisdescr,"...",sep=""))

      #---- Retrieve the variable in a scratch array. -------------------------------------#
      thisnow      = avenow[[thisvar]]

      #---- Find the variable range to make the colours look the same. --------------------#
      if (diffnow){
         thislevels   = pretty(c(-rangevar[vv,3],rangevar[vv,3]),n=ncolours)
         thisncolours = length(thislevels)
         thisrainbow  = thiscolour(n=thisncolours)
         #---- Make the values too close to zero to be white. -----------------------------#
         sel              = thislevels <   thissmfrac * rangevar[vv,3] &
                            thislevels > - thissmfrac * rangevar[vv,3]
         thisrainbow[sel] = "white"
      }else{
         thislevels   = pretty(c(rangevar[vv,1],rangevar[vv,2]),n=ncolours)
         thisncolours = length(thislevels)
         thisrainbow  = thiscolour(n=thisncolours)
      }#end if

      #----- Check whether the output directory exists. -----------------------------------#
      varpath=paste(figpath,thisoutpref,sep="/")
      if (! file.exists(varpath)) dir.create(varpath)


      #----- Find the title. --------------------------------------------------------------#
      wlaba = chron(round(data[[n]]$gtime[             1]),out.format=c(dates="mon day"))
      wlabz = chron(round(data[[n]]$gtime[data[[n]]$tmax]),out.format=c(dates="mon day"))
      letitre = paste(titlenow," - ",thisdescr,"\n Period: ",wlaba," - ",wlabz,", "
                     ,seasonnow,sep="")

      #------------------------------------------------------------------------------------#
      # Loop over the output formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in 1:noutput){

         #----- Find the output file. -----------------------------------------------------#
         outfile = paste(varpath,"/",thisoutpref,"-seasonmean-",suffoutnow,"."
                        ,plotform[o],sep="")

         #----- Open the appropriate device. ----------------------------------------------#
         if (plotform[o] == "x11"){
            X11(width=size$width,height=size$height,pointsize=ptsz)
         }else if(plotform[o] == "eps"){
             postscript(file=outfile,width=size$width,height=size$height,pointsize=ptsz,
                        paper=paper)
         }else{
             png(filename=outfile,width=size$width*depth,height=size$height*depth,
                 res=depth,pointsize=ptsz)
         }#end if

         #---------------------------------------------------------------------------------#
         #    Plot the variable.                                                           #
         #---------------------------------------------------------------------------------#
         sombreado(x=avenow$glon,y=avenow$glat,z=thisnow
                  ,levels=thislevels,nlevels=thisncolours,col=thisrainbow
                  ,plot.title=title(main=letitre,xlab="",ylab=""
                                   ,cex.main=cexaxes,font.sub=2)
                  ,key.title=title(main=thisunit,cex.main=1.0)
                  ,key.axes =axis(side=4,cex.axis=cexaxes)
                  ,plot.axes={ axis(side=1,cex.axis=cexaxes)
                               axis(side=2,cex.axis=cexaxes)
                               if (addgrid) grid(lty="dotted",col="lightgray")
                               if (plotmap){
                                  if (southamerica){
                                     southammap(mapdir=mapdir,col="black",lwd=2
                                               ,lty="dotted")
                                  }else{
                                     map("state",fill=FALSE,boundary=TRUE,add=TRUE
                                        ,col="gray",lty="dotted",lwd=2)
                                     map("worldHires",fill=FALSE,boundary=TRUE
                                        ,add=TRUE,col="black",lty="dotted",lwd=2)
                                  }#end if
                               }else if (plotrivernow[aa]){
                                  contour(x=avenow$glon,y=avenow$glat
                                         ,z=data[[n]]$land[1,1,,]
                                         ,levels=landmin,drawlabels=FALSE
                                         ,col="black",lty="dotted",add=TRUE,lwd=2)
                               } #end if

                               #----- Add the contours to tell where the forest is. -------#
                               if (plotdefnow[aa]){
                                  contour(x=avenow$glon,y=avenow$glat
                                         ,z=data[[aa]]$vegz0.ps[1,1,,]
                                         ,levels=defthreshold,drawlabels=FALSE
                                         ,col="gray23",lty="solid",add=TRUE,lwd=2)
                               }#end if

                               if (plotroad){
                                  roadmap(rpath=rpath,lty="solid",lwd=2,col="black")
                               }#end if

                               #----- Add some location info. -----------------------------#
                               for (ss in 1:nsites){
                                  thissite = site[[ss]]
                                  text(x=thissite$lon,y=thissite$lat,labels=thissite$iata
                                      ,cex=1.0,font=2,col="gray16",adj=c(0.5,0.5),srt=0)
                               }#end for
                             } #end plot axes
                   ) #end sombreado

         #----- Close device. -------------------------------------------------------------#
         if (plotform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
      }#end for
   }#end for
}#end for
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Create the plots of seasonal mean diurnal cycle for the three averages, and all       #
# variables.                                                                               #
#------------------------------------------------------------------------------------------#
for (aa in 1:nlist){

   print (paste(" + Plotting the seasonal mean diurnal cycle for ",seasonwhen[aa],"...",
                sep=""))

   #----- Retrieve information about this average. ----------------------------------------#
   avenow     = dmean[[aa]]
   titlenow   = titlelist[aa]
   plotlainow = plotlai[aa]
   diffnow    = difference[aa]
   suffoutnow = suffout[aa]
   seasonnow  = seasonwhen[aa]

   for (vv in 1:nvars){

      #---- Find the output-related variable information. ---------------------------------#
      thisvar     = vardf$name        [vv]
      thisoutpref = vardf$outpref     [vv]
      thisdescr   = vardf$description [vv]
      thisunit    = vardf$unit        [vv]
      thiscolour  = get(vardf$colour  [vv])
      thissmfrac  = as.numeric(vardf$smallfrac   [vv])


      print (paste("   - Variable: ",thisdescr,"...",sep=""))

      #---- Retrieve the variable in a scratch array. -------------------------------------#
      thisnow      = avenow[[thisvar]]

      #---- Find the variable range to make the colours look the same. --------------------#
      if (diffnow){
         thislevels   = pretty(c(-rangevar[vv,6],rangevar[vv,6]),n=ncolours)
         thisncolours = length(thislevels)
         thisrainbow  = thiscolour(n=thisncolours)
         #---- Make the values too close to zero to be white. -----------------------------#
         sel              = thislevels <   thissmfrac * rangevar[vv,6] &
                            thislevels > - thissmfrac * rangevar[vv,6]
         thisrainbow[sel] = "white"
      }else{
         thislevels   = pretty(c(rangevar[vv,4],rangevar[vv,5]),n=ncolours)
         thisncolours = length(thislevels)
         thisrainbow  = thiscolour(n=thisncolours)
      }#end if

      #----- Check whether the output directory exists. -----------------------------------#
      varpath=paste(figpath,thisoutpref,sep="/")
      if (! file.exists(varpath)) dir.create(varpath)

      #------------------------------------------------------------------------------------#
      #     Loop over all hours.                                                           #
      #------------------------------------------------------------------------------------#
      for (h in 1:length(avenow$hour)){
         #----- Find the hour. ------------------------------------------------------------#
         hh = avenow$hour[h]
         chh=substring(100+hh,2,3)


         #----- Find the title. -----------------------------------------------------------#
         wlaba = chron(round(data[[n]]$gtime[   1]),out.format=c(dates="mon day"))
         wlabz = chron(round(data[[n]]$gtime[data[[n]]$tmax]),out.format=c(dates="mon day"))
         letitre = paste(titlenow," - ",thisdescr,"\n Period: ",wlaba," - ",wlabz,", "
                        ,seasonnow," - ",chh," UTC",sep="")

         #---------------------------------------------------------------------------------#
         # Loop over the output formats.                                                   #
         #---------------------------------------------------------------------------------#
         for (o in 1:noutput){

            #----- Find the output file. --------------------------------------------------#
            outfile = paste(varpath,"/",thisoutpref,"-dcycle-",chh,"-",suffoutnow,"."
                           ,plotform[o],sep="")

            #----- Open the appropriate device. -------------------------------------------#
            if (plotform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(plotform[o] == "eps"){
                postscript(file=outfile,width=size$width,height=size$height,pointsize=ptsz,
                           paper=paper)
            }else{
                png(filename=outfile,width=size$width*depth,height=size$height*depth,
                    res=depth,pointsize=ptsz)
            }#end if

            #------------------------------------------------------------------------------#
            #    Plot the variable.                                                        #
            #------------------------------------------------------------------------------#
            sombreado(x=avenow$glon,y=avenow$glat,z=thisnow[h,,]
                     ,levels=thislevels,nlevels=thisncolours,col=thisrainbow
                     ,plot.title=title(main=letitre,xlab="",ylab=""
                                      ,cex.main=cexaxes,font.sub=2)
                     ,key.title=title(main=thisunit,cex.main=1.0)
                     ,key.axes =axis(side=4,cex.axis=cexaxes)
                     ,plot.axes={ axis(side=1,cex.axis=cexaxes)
                                  axis(side=2,cex.axis=cexaxes)
                                  if (addgrid) grid(lty="dotted",col="lightgray")
                                  if (plotmap){
                                     if (southamerica){
                                        southammap(mapdir=mapdir,col="black",lwd=2
                                                  ,lty="dotted")
                                     }else{
                                        map("state",fill=FALSE,boundary=TRUE,add=TRUE
                                           ,col="gray",lty="dotted",lwd=2)
                                        map("worldHires",fill=FALSE,boundary=TRUE
                                           ,add=TRUE,col="black",lty="dotted",lwd=2)
                                     }#end if
                                  }else if (plotrivernow[aa]){
                                     contour(x=avenow$glon,y=avenow$glat
                                            ,z=data[[n]]$land[1,1,,]
                                            ,levels=landmin,drawlabels=FALSE
                                            ,col="black",lty="dotted",add=TRUE,lwd=2)
                                  } #end if

                                  #----- Add the contours to tell where the forest is. ----#
                                  if (plotdefnow[aa]){
                                     contour(x=avenow$glon,y=avenow$glat
                                            ,z=data[[aa]]$vegz0.ps[1,1,,]
                                            ,levels=defthreshold,drawlabels=FALSE
                                            ,col="gray23",lty="solid",add=TRUE,lwd=2)
                                  }#end if

                                  if (plotroad){
                                     roadmap(rpath=rpath,lty="solid",lwd=2,col="black")
                                  }#end if

                                  #----- Add some location info. --------------------------#

                                  for (ss in 1:nsites){
                                     thissite = site[[ss]]
                                     text(x=thissite$lon,y=thissite$lat,labels=thissite$iata
                                         ,cex=1.0,font=2,col="gray16",adj=c(0.5,0.5),srt=0)
                                  }#end for
                                } #end plot axes
                      ) #end sombreado

            #----- Close device. ----------------------------------------------------------#
            if (plotform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         }#end for
      }#end for
   }#end for
}#end for
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

print(" + Plotting the diurnal cycle of the POIs")
for (ipy in 1:npois){
   pp = poi[[ipy]]
   print(paste("   - ",pp$lieu,"...",sep=""))

   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   for (hh in 1:nhov){

      #----- Retrieve variable information from the list. ---------------------------------#
      hovdinow    = hovdi[[hh]]
      vnames      = hovdinow$vnam  
      description = hovdinow$desc  
      lcolours    = hovdinow$colour
      prefix      = hovdinow$prefix
      theme       = hovdinow$theme 
      unit        = hovdinow$unit  
      legpos      = hovdinow$legpos
      plotit      = hovdinow$plt   
 
      if (plotit){


         #----- Define the number of layers. ----------------------------------------------#
         nlayers = length(vnames)

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir = paste(pp$outpref,"tsvar",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",theme,"Time series ..."))
         for (aa in 1:nlist){
            pdm     = pp$dmean[[aa]]
            ylimit  = range(pdm[vnames],na.rm=TRUE)

            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:noutput){
               fichier = paste(outdir,"/",prefix,"-",suffout[aa],"-",pp$suffix,"."
                              ,plotform[o],sep="")
               if(plotform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(plotform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(plotform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if

               letitre = paste(seasonwhen[aa]," Simulation - ",pp$lieu,
                               " \n"," Time series of mean diurnal cycle: ",theme,sep="")

               plot(x=pp$hour,y=pdm[[vnames[1]]],type="n",main=letitre,xlab="Hour (GMT)",
                    ylim=ylimit,ylab=paste("[",unit,"]",sep=""))
               if (hovgrid) grid(col="lightgray",lty="dotted")
               for (l in 1:nlayers){
                  points(x=pp$hour,y=pdm[[vnames[l]]],pch=16,col=lcolours[l]
                        ,cex=1.2,type="o",lwd=2)
               }#end for
               legend(x=legpos,inset=0.05,legend=description,col=lcolours,pch=16,lwd=2)
               if (plotform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for plotform
         }#end for, nlist
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#
