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
ctlfiles = c("abarca-_g1.ctl","bbarca-_g1.ctl")
season   = c("BARCA A","BARCA B")
slabel   = c("abarca","bbarca")

#----- List of variables to be plotted. ---------------------------------------------------#
vardf        = data.frame(t(rep(NA,times=5)))
names(vardf) = c("name","outpref","description","unit","colour")
vardf[ 1,] = c("hflxca"  ,"sensible","Sensible heat flux"           ,    "W/m2" ,"muitas" )
vardf[ 2,] = c("qwflxca" ,"latent"  ,"Evapotranspiration"           ,  "mm/day" ,"muitas" )
vardf[ 3,] = c("gpp"     ,"gpp"     ,"Gross Primary Production"     ,"kgC/m2/yr","iatlas" )
vardf[ 4,] = c("plresp"  ,"plresp"  ,"Plant Respiration"            ,"kgC/m2/yr","iatlas" )
vardf[ 5,] = c("resphet" ,"hetresp" ,"Heterotrophic Respiration"    ,"kgC/m2/yr","iatlas" )
vardf[ 6,] = c("nep"     ,"nep"     ,"Net Ecosystem Production"     ,"kgC/m2/yr","iatlas" )
vardf[ 7,] = c("rsnet"   ,"rsnet"   ,"Net shortwave Radiation"      ,     "W/m2","icloudy")
vardf[ 8,] = c("rlong"   ,"rlong"   ,"Net longwave Radiation"       ,     "W/m2","icloudy")
vardf[ 9,] = c("tcan.ps" ,"tcan"    ,"Canopy Air Temperature"       ,    " degC","muitas" )
vardf[10,] = c("rvcan.ps","rvcan"   ,"Canopy Air Vapor Mixing Ratio",    " g/kg","imuitas")
vardf[11,] = c("ecoresp" ,"ecoresp" ,"Ecosystem Respiration"        ,"kgC/m2/yr","iatlas" )
vardf[12,] = c("prate"   ,"prate"   ,"Precipitation rate"           , "kg/m2/hr","imuitas")
#------------------------------------------------------------------------------------------#



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
#    Define the sub-domain you want to use.  This is probably not going to be the exact    #
# domain, but it will be the closest it can be.                                            #
#------------------------------------------------------------------------------------------#
lona    = -78.9
lonz    = -35.1
lata    = -14.2
latz    =   9.4

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

addgrid        = FALSE     # Should I include a grid?
southamerica   = TRUE      # Is this a South America run (for map plotting)?
plotroad       = FALSE     # Should I include the roads?
ncolours       = 200       # Approximate number of colours for image plot.
colscheme      = "imuitas" # Colour scheme for precipitation rate. 
logpal         = FALSE     # Should the colour palette be in log scale? 
fontsub        = 2         # Type of font for subtitle:

smallfrac      = 0.05      # A small fraction of the range 
cexaxes        = 1.5       # Expansion factor for axes.

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
nvars = length(vardf[,1])

#----- List of sites. ---------------------------------------------------------------------#
nsites = length(site)

#----- List of CTL files. -----------------------------------------------------------------#
nctl   = length(ctlfiles)

#----- Load functions in Rsc. -------------------------------------------------------------#
source(paste(srcdir,"atlas.r",sep="/"))
source(paste(srcdir,"bintemp.r",sep="/"))
source(paste(srcdir,"cloudy.r",sep="/"))
source(paste(srcdir,"eps.r",sep="/"))
source(paste(srcdir,"gridp.r",sep="/"))
source(paste(srcdir,"gridt.r",sep="/"))
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

#----- Getting rid of the annoying bell. --------------------------------------------------#
options(locatorBell=FALSE)

#----- Make plotform case insensitive. ----------------------------------------------------#
plotform = tolower(plotform)

#----- Determine the number of plot output formats. ---------------------------------------#
noutput = length(plotform)

#----- Find the time for plot. ------------------------------------------------------------#
whena = chron(dates=paste(montha,daya,yeara,sep="/"),times=paste(houra,0,0,sep=":"))
whenz = chron(dates=paste(monthz,dayz,yearz,sep="/"),times=paste(hourz,0,0,sep=":"))

#----- The list of variables to be loaded from the file. ----------------------------------#
vlist = c("lai.ps","agb.ps","vegz0.ps","hflxca","qwflxca","rshort","albedt","rlong","gpp"
         ,"plresp","resphet","tcan.ps","rvcan.ps","conprr01","pcprate")

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
   namesaver  = list()
   dmean      = list()
   namesdmean = list()

   #----- Read the data from all simulations. ---------------------------------------------#
   print (" + Reading BRAMS grid...")
   for (n in 1:nctl){
      data[[n]] = readgrads(vari=vlist,info=info[[n]],xlim=c(xa,xz),ylim=c(ya,yz),
                          tlim=c(ta[n],tz[n]))
      #------------------------------------------------------------------------------------#
      #     Compute/rescale additional variables.                                          #
      #------------------------------------------------------------------------------------#
      #----- Net shortwave radiation. -----------------------------------------------------#
      data[[n]]$rsnet = data[[n]]$rshort*(1.-data[[n]]$albedt)
      #----- Convert latent heat to evaporation in mm/day. --------------------------------#
      data[[n]]$qwflxca = data[[n]]$qwflxca * day.sec / alvl
      #----- Convert the ecosystem properties to kgC/m2/yr. -------------------------------#
      data[[n]]$gpp     = data[[n]]$gpp     * umols.2.kgCyr
      data[[n]]$plresp  = data[[n]]$plresp  * umols.2.kgCyr
      data[[n]]$resphet = data[[n]]$resphet * umols.2.kgCyr
      #----- Find the ecosystem respiration. ----------------------------------------------#
      data[[n]]$ecoresp = data[[n]]$plresp + data[[n]]$resphet
      #----- Find NEP. --------------------------------------------------------------------#
      data[[n]]$nep     = data[[n]]$gpp - data[[n]]$plresp - data[[n]]$resphet
      #----- Find total precipitation rate. -----------------------------------------------#
      data[[n]]$prate   = data[[n]]$conprr01 + data[[n]]$pcprate
      #------------------------------------------------------------------------------------#

      #----- Find the hours. --------------------------------------------------------------#
      data[[n]]$hour  = hours(data[[n]]$gtime)
      #------------------------------------------------------------------------------------#



      #----- Read the data from all simulations. ------------------------------------------#
      allhrs          = sort(unique(data[[n]]$hour))
      aver[[n]]       = list(xmax=data[[n]]$xmax,ymax=data[[n]]$ymax,zmax=data[[n]]$zmax
                            ,glev=data[[n]]$glev,glon=data[[n]]$glon,glat=data[[n]]$glat)
      namesaver[[n]]  = names(aver[[n]])
      dmean[[n]]      = list(xmax=data[[n]]$xmax,ymax=data[[n]]$ymax,zmax=data[[n]]$zmax
                            ,tmax=length(allhrs),hour=allhrs,glev=data[[n]]$glev
                            ,glon=data[[n]]$glon,glat=data[[n]]$glat)
      namesdmean[[n]] = names(dmean[[n]])
      #------------------------------------------------------------------------------------#
   }#end for
   aver[[3]]       = aver[[1]]
   dmean[[3]]      = dmean[[1]]
   namesaver[[3]]  = namesaver[[1]]
   namesdmean[[3]] = namesdmean[[1]]

   #----- Loop over all variables to compute the diurnal cycle. ---------------------------#
   for (vv in 1:nvars){
      thisvar   = vardf$name[vv]
      thisdescr = vardf$description[vv]
      print (paste(" + Computing the seasonal average of ",thisdescr,"...",sep=""))

      #----- Copy the variable to a scratch matrix. ---------------------------------------#
      aux      = list()
      for (n in 1:2){
         aux[[n]] = matrix(data[[n]][[thisvar]],nrow=data[[n]]$tmax,
                                                ncol=data[[n]]$xmax*data[[n]]$ymax)
         seasmean = apply(X=aux[[n]]  ,MARGIN=2,FUN=mean,na.rm=TRUE)
         dcycmean = qapply(mat=aux[[n]],bycol=TRUE,index=data[[n]]$hour
                          ,func=mean,na.rm=TRUE)

         #----- Compute the mean differences, then bind the differences to the list. ------#
         aver[[n]]$this  = matrix(seasmean,nrow=aver[[n]]$xmax  ,ncol=aver[[n]]$ymax  )
         dmean[[n]]$this = array(dcycmean,
                                 dim=c(dmean[[n]]$tmax,dmean[[n]]$xmax,dmean[[n]]$ymax))
      }#end for
      aver[[3]]$this  = aver[[2]]$this  - aver[[1]]$this
      dmean[[3]]$this = dmean[[2]]$this - dmean[[1]]$this

      for (n in 1:3){
         #----- Update the list names. ----------------------------------------------------#
         namesaver[[n]]    = c(namesaver[[n]]  ,thisvar)
         names(aver[[n]])  = namesaver[[n]]

         namesdmean[[n]]   = c(namesdmean[[n]] ,thisvar)
         names(dmean[[n]]) = namesdmean[[n]]
      }#end for
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

#----- Define some instructions on what and how to plot the averages. ---------------------#
nave       = length(aver)
titlelist  = c("Seasonal mean","Seasonal mean","Mean difference")
plotlai    = c(FALSE, TRUE,TRUE)
difference = c(FALSE,FALSE,TRUE)
seasonwhen = c(season,paste(season[2],season[1],sep=" - "))
suffout    = c(slabel,"diff")

#------------------------------------------------------------------------------------------#
#    Create the output for the three averages, and all variables.                          #
#------------------------------------------------------------------------------------------#
for (aa in 1:nave){

   print (paste(" + Plotting averages for ",seasonwhen[aa],"...",sep=""))

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


      print (paste("   - Plotting the seasonal mean of ",thisdescr,"...",sep=""))

      #---- Retrieve the variable in a scratch array. -------------------------------------#
      thisnow      = avenow[[thisvar]]

      #---- Find the variable range to make the colours look the same. --------------------#
      if (diffnow){
         thislevels   = pretty(c(-rangevar[vv,3],rangevar[vv,3]),n=ncolours)
         thisncolours = length(thislevels)
         thisrainbow  = thiscolour(n=thisncolours)
         #---- Make the values too close to zero to be white. -----------------------------#
         sel              = thislevels <   smallfrac * rangevar[vv,3] &
                            thislevels > - smallfrac * rangevar[vv,3]
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
      letitre=paste(titlenow," - ",thisdescr,"\n Period: ",seasonnow,sep="")

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
                               if (southamerica){
                                  southammap(mapdir=mapdir,col="black")
                               }else{
                                  map("state",fill=FALSE,boundary=TRUE,add=TRUE
                                     ,col="gray")
                                  map("worldHires",fill=FALSE,boundary=TRUE,add=TRUE
                                     ,col="Black")
                               } #end if

                               if (plotroad){
                                  roadmap(rpath=rpath,lty="solid",lwd=2,col="yellow2")
                               }#end if

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
#    Create the output for the three averages, and all variables.                          #
#------------------------------------------------------------------------------------------#
for (aa in 1:nave){

   print (paste(" + Plotting averages for ",seasonwhen[aa],"...",sep=""))

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


      print (paste("   - Plotting the mean diurnal cycle of ",thisdescr,"...",sep=""))

      #---- Retrieve the variable in a scratch array. -------------------------------------#
      thisnow      = avenow[[thisvar]]

      #---- Find the variable range to make the colours look the same. --------------------#
      if (diffnow){
         thislevels   = pretty(c(-rangevar[vv,6],rangevar[vv,6]),n=ncolours)
         thisncolours = length(thislevels)
         thisrainbow  = thiscolour(n=thisncolours)
         #---- Make the values too close to zero to be white. -----------------------------#
         sel              = thislevels <   smallfrac * rangevar[vv,6] &
                            thislevels > - smallfrac * rangevar[vv,6]
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
         print(paste("     * ",chh," UTC",sep=""))

         #----- Find the title. -----------------------------------------------------------#
         letitre = paste(titlenow," - ",thisdescr,"\n Period: ",seasonnow," - ",chh
                        ," UTC",sep="")
      

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
                                  if (southamerica){
                                     southammap(mapdir=mapdir,col="black")
                                  }else{
                                     map("state",fill=FALSE,boundary=TRUE,add=TRUE
                                        ,col="gray")
                                     map("worldHires",fill=FALSE,boundary=TRUE,add=TRUE
                                        ,col="Black")
                                  } #end if

                                  if (plotroad){
                                     roadmap(rpath=rpath,lty="solid",lwd=2,col="yellow2")
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
