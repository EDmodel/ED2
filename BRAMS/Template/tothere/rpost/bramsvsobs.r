#------------------------------------------------------------------------------------------#
#     Leave the following lines as the first commands.  DO NOT define anything before this #
# line.                                                                                    #
#------------------------------------------------------------------------------------------#
graphics.off()
readitagain = ! "brams" %in% ls()
if (readitagain){
   rm(list=ls())
   readdata = TRUE
}else{
   readdata = FALSE
}#end if

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
ctlfile = "owetamc-_g2.ctl"                # WETAMC file
prefout = "wetamc"                         # Prefix for output figures
dataset = "BRAMS"                          # Dataset used in the plots (for title in
                                           #    figures only)
period  = "WETAMC"                         # Period (for title in figures only)
trmmnw  = 9                                # Number of this network. Any number is fine...
iata    = "BRM"                            # 3-letter identifier for this data.
dtime   = 1                                # Nominal delta-t (in hours)
path    = getwd()                          # Working directory.
libcomp = paste(getwd(),"libcomp",sep="/")


#------------------------------------------------------------------------------------------#
#    Define the sub-domain you want to use.  This is probably not going to be the exact    #
# domain, but it will be the closest it can be.                                            #
#------------------------------------------------------------------------------------------#
lona    = -67.0
lonz    = -57.0
lata    = -15.0
latz    =  -7.0

#----- Initial time for diurnal property analysis. ----------------------------------------#
yeara   =  1999        # Year
montha  =     1        # Month
daya    =    20        # Day
houra   =     0        # Hour

#----- Final time for diurnal property analysis. ------------------------------------------#
yearz   =  1999        # Year
monthz  =     2        # Month
dayz    =    20        # Day
hourz   =     0        # Hour
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------# 
#    Define the "chron" origin.  This only matters if your dataset has days in day of year #
# (doy, or often incorrectly called "Julian day").                                         #
#    In R, Julian day is the day "0", so if your day 1 means January 1st, 1999, for        #
# instance, then set your chron origin to December 31, 1998.                               #
#------------------------------------------------------------------------------------------#
yearo  = 1998
montho =   12
dayo   =   31

#------------------------------------------------------------------------------------------#
#    Define the paths for all datasets we will look at.                                    #
# trmmpath  - Directory where the TRMM statistics files are.                               #
# rgpath    - Directory where the rain gauge statistics files are.                         #
# awspath   - Directory where the automatic weather station statistics files are.          #
# radiopath - Directory where the radiosond 
#------------------------------------------------------------------------------------------#
trmmpath  = "/n/moorcroft_scratch/mlongo/lba_data/wetamc/TRMM/output/"
rgpath    = "/n/moorcroft_scratch/mlongo/lba_data/wetamc/rain_gauge/output/"
awspath   = "/n/moorcroft_scratch/mlongo/lba_data/wetamc/aws/stats/"
radiopath = "/n/moorcroft_scratch/mlongo/lba_data/wetamc/radio/stats/"

#----- NA for radiosondes. ----------------------------------------------------------------#
na.radio = c("NA","-9999.99")

#----- Output options. --------------------------------------------------------------------#
plotform = c("png")#,"x11")  # Type of output: currently valid options are png, eps, and x11.
paper    = "letter"        # Paper size
depth    =   96            # Resolution in pixels per inch (png only)
ptsz     =   12            # Font size
offnet   = 0.04            # Offset for network precipitation map.
offdcyc  = 0.60            # Offset for the diurnal cycle plots.

addgrid        = TRUE      # Should I include a grid?
southamerica   = TRUE      # Is this a South America run (for map plotting)?
ncols          = 200       # Approximate number of colours for image plot.
colscheme      = "imuitas" # Colour scheme for precipitation rate. 
logpal         = FALSE     # Should the colour palette be in log scale? 
fontsub        = 2         # Type of font for subtitle:
                           # 1. Regular; 2. Bold; 3. Italic; 4. Bold-Italic.


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#






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
timeatstart=proc.time()[3]

print (" + Initial settings...")

#----- Load libraries. --------------------------------------------------------------------#
isok = require(chron)
isok = require(fields)
isok = require(akima)

#----- Define useful paths. ---------------------------------------------------------------#
mapdir = paste(path,"samap",sep="/")
srcdir = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc"

#----- Getting rid of the annoying bell. --------------------------------------------------#
options(locatorBell=FALSE)

#----- Make plotform case insensitive. ----------------------------------------------------#
plotform = tolower(plotform)

#----- Determine the number of plot output formats. ---------------------------------------#
noutput = length(plotform)

#----- Find the initial and final time for analysis. --------------------------------------#
whena = chron(dates=paste(montha,daya,yeara,sep="/"),times=paste(houra,0,0,sep=":"))
whenz = chron(dates=paste(monthz,dayz,yearz,sep="/"),times=paste(hourz,0,0,sep=":"))

#----- Load functions in Rsc. -------------------------------------------------------------#
source(paste(srcdir,"bintemp.r",sep="/"))
source(paste(srcdir,"bramsvsaws.r",sep="/"))
source(paste(srcdir,"bramsvsradio.r",sep="/"))
source(paste(srcdir,"cut2d.r",sep="/"))
source(paste(srcdir,"bramsvsrg.r",sep="/"))
source(paste(srcdir,"bramsvstrmm.r",sep="/"))
source(paste(srcdir,"qapply.r",sep="/"))
source(paste(srcdir,"eps.r",sep="/"))
source(paste(srcdir,"gridp.r",sep="/"))
source(paste(srcdir,"gridt.r",sep="/"))
source(paste(srcdir,"interpol.r",sep="/"))
source(paste(srcdir,"lnlikeaws.r",sep="/"))
source(paste(srcdir,"lnlikeradio.r",sep="/"))
source(paste(srcdir,"lnlikerg.r",sep="/"))
source(paste(srcdir,"lnliketrmm.r",sep="/"))
source(paste(srcdir,"muitas.r",sep="/"))
source(paste(srcdir,"nearest.neighbour.r",sep="/"))
source(paste(srcdir,"plotsize.r",sep="/"))
source(paste(srcdir,"pretty.log.r",sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
source(paste(srcdir,"readctl.r",sep="/"))
source(paste(srcdir,"readgrads.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"thermlib.r",sep="/"))


#----- Plot size. -------------------------------------------------------------------------#
size = plotsize(proje=TRUE,limlon=c(lona,lonz),limlat=c(lata,latz),extendfc=TRUE,
                paper=paper)


#------------------------------------------------------------------------------------------#
#      Determine the output and plot output directories, and check whether they exist or   #
# not.  In case they don't, create them.                                                   #                
#------------------------------------------------------------------------------------------#
print (" + Checking output paths...")
figpath = paste(path,"figures",sep="/")
outpath = paste(path,"output",sep="/")
if (! file.exists(figpath)) dir.create(figpath)
if (! file.exists(outpath)) dir.create(outpath)


#----- This part is done only when reading the data is required. --------------------------#
if (readdata){
   #---------------------------------------------------------------------------------------#
   #     Retrieve the summaries of all datasets.  They contain the coordinates and point   #
   # names, among other things...                                                          #
   #---------------------------------------------------------------------------------------#
   print (" + Reading the observation summaries...")
   #----- TRMM. ---------------------------------------------------------------------------#
   print ("   - TRMM summary...")
   trmmsmry        = read.table(file=paste(trmmpath,"trmmnet-smry.txt",sep="/"),
                                col.names=c("wlon","elon","slat","nlat","net","accr"))
   trmmsmry$swedge = cbind(trmmsmry$wlon,trmmsmry$slat)
   trmmsmry$needge = cbind(trmmsmry$elon,trmmsmry$nlat)
   #----- Rain gauge network. -------------------------------------------------------------#
   print ("   - Rain gauge network summary...")
   rgsmry      = read.table(file=paste(rgpath,"rgnet-smry.txt",sep="/"),
                            col.names=c("lon","lat","net","name","accr"))
   #----- AWS network. --------------------------------------------------------------------#
   print ("   - Automatic weather station network summary...")
   awssmry   = read.table(file=paste(awspath,"aws-smry.txt",sep="/"),
                          col.names=c("lon","lat","name"))
   #----- Radiosonde network. -------------------------------------------------------------#
   print ("   - Radiosonde network summary...")
   radiosmry        = read.table(file=paste(radiopath,"radio-smry.txt",sep="/"),
                                 col.names=c("lon","lat","name"))
   radiosmry$nradio = dim(radiosmry)[1]

   obs = list()
   #----- Read TRMM statistics and store them on a list named TRMM. -----------------------#
   print (" + Reading TRMM statistics...")
   obs$trmm = list()
   obs$trmm$dmean  = read.table(file=paste(trmmpath,"trmmnet-mean.txt",sep="/"),header=TRUE)
   obs$trmm$dsdev  = read.table(file=paste(trmmpath,"trmmnet-sdev.txt",sep="/"),header=TRUE)
   obs$trmm$dprob  = read.table(file=paste(trmmpath,"trmmnet-prob.txt",sep="/"),header=TRUE)
   obs$trmm$dccnt  = read.table(file=paste(trmmpath,"trmmnet-ccnt.txt",sep="/"),header=TRUE)
   obs$trmm$xhday = names(obs$trmm$dmean)

   #----- Read the rain gauge statistics and store them on a list named rg. ---------------#
   print (" + Reading rain gauge statistics...")
   obs$rg = list()
   obs$rg$dmean  = read.table(file=paste(rgpath,"rgnet-mean.txt",sep="/"),header=TRUE)
   obs$rg$dsdev  = read.table(file=paste(rgpath,"rgnet-sdev.txt",sep="/"),header=TRUE)
   obs$rg$dprob  = read.table(file=paste(rgpath,"rgnet-rfrq.txt",sep="/"),header=TRUE)
   obs$rg$dccnt  = read.table(file=paste(rgpath,"rgnet-ccnt.txt",sep="/"),header=TRUE)
   obs$rg$xhday  = names(obs$rg$dmean)
   obs$rg$nrgs   = dim(rgsmry)[1]

   #----- Read automatic weather station statistics and store them on a list named aws. ---#
   print (" + Reading AWS statistics...")
   obs$aws = list()
   obs$aws$naws  = dim(awssmry)[1]
   obs$aws$vars = c("psfc","thsfc","rvsfc","wsfc")
   obs$aws$nvars = length(obs$aws$vars)
   obs$aws$places = awssmry$name
   for(vars in obs$aws$vars){
      obs$aws[[vars]] = list()
   }

   first = TRUE
   for (a in 1: obs$aws$naws){

      #----- Read the statistics. ---------------------------------------------------------#
      awsdmean     = paste(awspath,"/",awssmry$name[a],"-aws-mean.txt",sep="")
      dmean        = read.table(file=awsdmean,header=TRUE,row.names=1)
      names(dmean) = obs$aws$vars
      awsdsdev     = paste(awspath,"/",awssmry$name[a],"-aws-sdev.txt",sep="")
      dsdev        = read.table(file=awsdsdev,header=TRUE,row.names=1)
      names(dsdev) = obs$aws$vars
      awsdccnt     = paste(awspath,"/",awssmry$name[a],"-aws-ccnt.txt",sep="")
      dccnt        = read.table(file=awsdccnt,header=TRUE,row.names=1)
      names(dccnt) = obs$aws$vars

      if (first){
         #----- Allocate the common data frame. -------------------------------------------#
         obs$aws$xhday = row.names(dmean)
         if (substr(obs$aws$xhday[1],1,1) != "X"){
            obs$aws$xhday = paste("X",obs$aws$xhday,sep="")
         }#end if
         matini = matrix(NA,ncol=obs$aws$naws,nrow=length(obs$aws$xhday),
                            dimnames=list(obs$aws$xhday,obs$aws$places))
         for (vars in obs$aws$vars){
            obs$aws[[vars]]$dmean = matini
            obs$aws[[vars]]$dsdev = matini
            obs$aws[[vars]]$dccnt = matini
         }#end for
         rm(matini)
         first = FALSE
      }#end first

      #------------------------------------------------------------------------------------#
      #      Copy variables to the structures.                                             #
      #------------------------------------------------------------------------------------#
      for (vari in obs$aws$vars){
         obs$aws[[vari]]$dmean[,a]  = dmean[,vari]
         obs$aws[[vari]]$dsdev[,a]  = dsdev[,vari]
         obs$aws[[vari]]$dccnt[,a]  = dccnt[,vari]
      }#end for
   }#end for

   #----- Read automatic weather station statistics and store them on a list named aws. ---#
   print (" + Reading the radiosonde statistics...")
   obs$radio = list()
   obs$radio$nradio  = dim(radiosmry)[1]
   obs$radio$vars    = c("zgeo","theta","rvap","uspd","vspd")
   obs$radio$nvars   = length(obs$radio$vars)
   obs$radio$places  = radiosmry$name
   for(vars in obs$radio$vars){
      obs$radio[[vars]] = list()
   }#end for (vars in obs$radio$vars)

   first = TRUE
   for (a in 1: obs$radio$nradio){
      iata = radiosmry$name[a]
      for (vari in obs$radio$vars){

         #----- Read the statistics. ------------------------------------------------------#
         radiopref = paste(radiopath,"/obs-",iata,"-",sep="")

         radiodmean   = paste(radiopref,"mean-",vari,".txt",sep="")
         dmean        = read.table(file=radiodmean,header=TRUE,row.names=1,
                                   na.strings=na.radio)
         radiodsdev   = paste(radiopref,"sdev-",vari,".txt",sep="")
         dsdev        = read.table(file=radiodsdev,header=TRUE,row.names=1,
                                   na.strings=na.radio)
         radiodccnt   = paste(radiopref,"ccnt-",vari,".txt",sep="")
         dccnt        = read.table(file=radiodccnt,header=TRUE,row.names=1,
                                   na.strings=na.radio)
         if (first){
            #----- Allocate the common data frame. ----------------------------------------#
            obs$radio$xhday  = row.names(dmean)
            obs$radio$plevs  = names(dmean)
            obs$radio$nplevs = length(obs$radio$plevs)
            if (substr(obs$radio$xhday[1],1,1) != "X"){
               obs$radio$xhday = paste("X",obs$radio$xhday,sep="")
            }#end if
            if (substr(obs$radio$plevs[1],1,1) == "X"){
               obs$radio$plevs = as.integer(substring(obs$radio$plevs,2))
            }
            arrini = array(NA,dim=c(obs$radio$nradio
                                   ,length(obs$radio$xhday)
                                   ,length(obs$radio$plevs)))
            for (varvar in obs$radio$vars){
               dimnames(arrini) = list(radiosmry$name,obs$radio$xhday,obs$radio$plevs)
               obs$radio[[varvar]]$dmean = arrini
               obs$radio[[varvar]]$dsdev = arrini
               obs$radio[[varvar]]$dccnt = arrini
            }#end for
            rm(arrini)
            first = FALSE
         }#end first
         #---------------------------------------------------------------------------------#
         #      Copy variables to the structures.                                          #
         #---------------------------------------------------------------------------------#

         obs$radio[[vari]]$dmean[iata,,]  = as.matrix(dmean)
         obs$radio[[vari]]$dsdev[iata,,]  = as.matrix(dsdev)
         obs$radio[[vari]]$dccnt[iata,,]  = as.matrix(dccnt)
      }#end for
   }#end for
   #----- Switch order of variable dimensions. --------------------------------------------#
   for (vari in obs$radio$vars){
      obs$radio[[vari]]$dmean = aperm(a=obs$radio[[vari]]$dmean,perm=c(2,3,1))
      obs$radio[[vari]]$dsdev = aperm(a=obs$radio[[vari]]$dsdev,perm=c(2,3,1))
      obs$radio[[vari]]$dccnt = aperm(a=obs$radio[[vari]]$dccnt,perm=c(2,3,1))
   }#end for

   #----- Retrieving the control file. ----------------------------------------------------#
   print (" + Reading BRAMS control file...")
   info = readctl(ctlfile)

   #----- Find the domain boundaries. -----------------------------------------------------#
   print (" + Finding the domain boundaries...")
   xa = which.min(abs(info$glon-lona))
   xz = which.min(abs(info$glon-lonz))
   ya = which.min(abs(info$glat-lata))
   yz = which.min(abs(info$glat-latz))
   ta = which.min(abs(info$gtime-whena))
   tz = which.min(abs(info$gtime-whenz))

   #----- Read the precipitation rate. ----------------------------------------------------#
   print (" + Reading BRAMS grid...")
   bramsgrid = readgrads(vari=c("lon","lat","topo"),info=info,xlim=c(xa,xz),ylim=c(ya,yz),
                         tlim=c(ta,ta))

   bramsgrid$point = data.frame(x=as.vector(bramsgrid$lon),y=as.vector(bramsgrid$lat))

   #----- Find the index of each TRMM dataset associated with each BRAMS grid point. ------#
   print (" + Mapping the BRAMS grid to the TRMM dataset...")
   xbreaks  = unique(sort(c(trmmsmry$swedge[,1],trmmsmry$needge[,1])))
   ybreaks  = unique(sort(c(trmmsmry$swedge[,2],trmmsmry$needge[,2])))
   vccut2d  = cut2d(x=bramsgrid$point,xbreaks=xbreaks,ybreaks=ybreaks)
   vclevs2d = levels(vccut2d)
   bramsgrid$maptrmm = match(vccut2d,vclevs2d)

   print (" + Finding the BRAMS points that are the closest to each rain gauge...")
   rgcoord = cbind(rgsmry$lon,rgsmry$lat)
   bramsgrid$nnrg = nearest.neighbour(netwcoord=rgcoord,gridcoord=bramsgrid$point)

   print (" + Finding the BRAMS points that are the closest to each AWS...")
   awscoord = cbind(awssmry$lon,awssmry$lat)
   bramsgrid$nnaws = nearest.neighbour(netwcoord=awscoord,gridcoord=bramsgrid$point)

   print (" + Finding the BRAMS points that are the closest to each radiosonde...")
   radiocoord = cbind(radiosmry$lon,radiosmry$lat)
   bramsgrid$nnradio = nearest.neighbour(netwcoord=radiocoord,gridcoord=bramsgrid$point)

   print (" + Reading BRAMS variables...")
   bramsraw = readgrads(vari=c("conprr01","pcprate","pcan.ps","theta2m","rv2m","u10m",
                               "press","theta","rv","ue.avg","ve.avg"),
                        info=info,xlim=c(xa,xz),ylim=c(ya,yz),tlim=c(ta,tz))

   #----- "Create" the geopotential variable. ---------------------------------------------#
   topo = as.vector(rep(bramsgrid$topo,times=info$lmax))
   rtgt = as.vector(rep(1. - bramsgrid$topo / info$glev[max(info$zmax)],times=info$lmax))
   sigz = as.vector(rep(info$glev,each=bramsgrid$xmax*bramsgrid$ymax))
   zgeo = rep(topo + sigz * rtgt,times=bramsgrid$tmax)
   zgeo = array(zgeo,dim=c(bramsgrid$xmax,bramsgrid$ymax,info$lmax,bramsraw$tmax))
   zgeo = aperm(a=zgeo,perm=c(4,3,1,2))
   bramsraw$zgeo = array(zgeo,dim=dim(bramsraw$press))

   #----- Change the name of surface variables so they have the same name as observations. #
   print (" + Standardising the surface variable names...")
   namesbramsraw = names(bramsraw)
   iswitch = match(c("pcan.ps","theta2m","rv2m","u10m"),namesbramsraw)
   namesbramsraw[iswitch] = obs$aws$vars
   names(bramsraw) = namesbramsraw

   #----- Change the name of level variables so they have the same name as observations. --#
   print (" + Standardising the height-level variable names...")
   namesbramsraw = names(bramsraw)
   iswitch = match(c("zgeo","theta","rv","ue.avg","ve.avg"),namesbramsraw)
   namesbramsraw[iswitch] = obs$radio$vars
   names(bramsraw) = namesbramsraw

   #---------------------------------------------------------------------------------------#
   #     BRAMS is a list that will contain the variables that matter for the support       #
   # calculation, so they will be pre-processed to have the same information as the        #
   # dataset to which the model will be compared.                                          #
   #---------------------------------------------------------------------------------------#
   brams = list()

   #---------------------------------------------------------------------------------------#
   #       Compute the precipitation rates.                                                #
   #---------------------------------------------------------------------------------------#
   print (" + Finding the total precipitation rate from BRAMS...")
   bramsraw$prate     = bramsraw$conprr01 + bramsraw$pcprate

   #----- 1. Comparison between BRAMS and TRMM. -------------------------------------------#
   print(" + Building the statistics to compare BRAMS and TRMM...")
   brams$trmm = bramsvstrmm(bramsgrid=bramsgrid,bramsraw=bramsraw,obs=obs)

   #----- 2. Comparison between BRAMS and rain gauge. -------------------------------------#
   print(" + Building the statistics to compare BRAMS and the rain gauge network...")
   brams$rg   = bramsvsrg(bramsgrid=bramsgrid,bramsraw=bramsraw,obs=obs)

   #----- 3. Comparison between BRAMS and AWS. --------------------------------------------#
   print(" + Building the statistics to compare BRAMS and the AWS network...")
   brams$aws  = bramsvsaws(bramsgrid=bramsgrid,bramsraw=bramsraw,obs=obs)

   #----- 4. Comparison between BRAMS and Radiosonde. -------------------------------------#
   print(" + Building the statistics to compare BRAMS and the radiosonde network...")
   brams$radio  = bramsvsradio(bramsgrid=bramsgrid,bramsraw=bramsraw,obs=obs)

}#end if readdata

#------------------------------------------------------------------------------------------#
#     At this point, brams and obs have structures with the same dimensions, so we can     #
# compute the support, which is the logarithm of likelihood.  support is a list with the   #
# likelihood of each dataset, computed for each time of the day.  The aggregated support   #
# will be computed later.                                                                  #
#------------------------------------------------------------------------------------------#
suppaux=list()

#----- 1. Compute the support between BRAMS and TRMM. -------------------------------------#
print(" + Computing the support between BRAMS and TRMM...")
suppaux$trmm  = lnliketrmm(brams=brams,obs=obs)

#----- 2. Compute the support between BRAMS and rain gauge. -------------------------------#
print(" + Computing the support between BRAMS and rain gauge network...")
suppaux$rg    = lnlikerg(brams=brams,obs=obs)

#----- 3. Compute the support between BRAMS and AWS. --------------------------------------#
print(" + Computing the support between BRAMS and automatic weather station network...")
suppaux$aws   = lnlikeaws(brams=brams,obs=obs)

#----- 4. Compute the support between BRAMS and Radiosonde. -------------------------------#
print(" + Computing the support between BRAMS and radiosonde network...")
suppaux$radio = lnlikeradio(brams=brams,obs=obs)

#------------------------------------------------------------------------------------------#
#    Compute the aggregated support, which encompasses all variables.  Each property       #
# (temperature, humidity, precipitation, etc.) for each data source (TRMM, radiosonde,     #
# etc.) will be treated as a dataset.  The total support is going to be the average of     #
# supports for each data source, which, in turn, is the average amongst all data.          #
#------------------------------------------------------------------------------------------#
print(" + Computing the aggregated support for each data source...")
aggr = list()
print("   - TRMM...")
aggr$trmm = mean(as.vector(as.matrix(suppaux$trmm)),na.rm=TRUE)
print("   - Rain gauge network...")
aggr$rg   = mean(as.vector(as.matrix(suppaux$rg))  ,na.rm=TRUE)
print("   - Automatic weather station...")
aggr$aws  = rep(NA,times=obs$aws$nvars)
for (v in 1:obs$aws$nvars){
   vari        = obs$aws$vars[v]
   aggr$aws[v] = mean(as.vector(suppaux$aws[[vari]]),na.rm=TRUE)
}#end for
print("   - Radiosonde...")
aggr$radio = rep(NA,times=obs$radio$nvars)
for (v in 1:obs$radio$nvars){
   vari          = obs$radio$vars[v]
   aggr$radio[v] = mean(as.vector(suppaux$radio[[vari]]),na.rm=TRUE)
}#end for
aggr = sapply(X=aggr,FUN=mean,na.rm=TRUE)

support = mean(aggr)

print (paste(" + The support for this data set is :",signif(support,5),"...",sep=""))

timeatend=proc.time()[3]
elapsedtime = timeatend - timeatstart

print(paste("=== Post-processing ends.  Run time: ",signif(elapsedtime,5),"... ===",sep=""))
