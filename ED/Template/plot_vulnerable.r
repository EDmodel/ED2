#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
here        = getwd()                          #   Current directory
srcdir      = "/n/home00/mlongo/util/Rsc"      #   Script directory
ibackground = 0                                #   Sought background (actual background
                                               #     is always transparent):
                                               #   0 -- White
                                               #   1 -- Black
                                               #   2 -- Dark grey
outroot     = file.path(here,paste("vulnerable_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Data set options.                                                                   #
#------------------------------------------------------------------------------------------#
retrieve.data    = TRUE                                            # Use previous Rdata
rdata.path       = file.path(here,"RData_vulnerable")              # Path with data.
rdata.vulnerable = file.path(rdata.path,"vulnerable_region.RData") # File name.
rdata.gridded    = file.path(rdata.path,"clim_region.RData")       # Gridded data set.
rdata.fini       = c(2008,12)                                      # Last time
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform         = c("png","pdf")  # Formats for output file.  Supported formats are:
                                  #   - "X11" - for printing on screen
                                  #   - "eps" - for postscript printing
                                  #   - "png" - for PNG printing
                                  #   - "pdf" - for PDF printing
depth           = 96              # PNG resolution, in pixels per inch
paper           = "letter"        # Paper size, to define the plot shape
wpaper          = "long"          # Paper size, to define the plot shape
ptsz            = 16              # Font size.
lwidth          = 2.5             # Line width
plotgrid        = TRUE            # Should I plot the grid in the background? 
n.colourbar     = 64              # Number of colours for the colour bars
theta           = 306.            # Azimuth for perspective projection
phi             = 25.             # Vertical angle for perspective projection
ltheta          = -210.           # Azimuth angle for light
shade           = 0.0001          # Shade intensity
expz            = 0.7             # Expansion factor for Z axis
proj.wall       = TRUE            # Project the points onto the wall
cex.pt3d        = 0.6             # Size for points in the 3-D plot
pch.skill       = 1               # Type of point to use in skill
pch.wmo         = 16              # Type of point to use in wmo plots
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Options for computing the vulnerability.                                             #
#------------------------------------------------------------------------------------------#
ntry.max        = 20              # Maximum number of attempts to find the score shift.
loud            = FALSE           # loud (verbose).
wmo.nyears.min  = 25              # Minimum number of years to consider the station valid
n.lsq.min       = 5               # Minimum number of years to consider in the LSQ
#------------------------------------------------------------------------------------------#




#----- Data sets. -------------------------------------------------------------------------#
datrain      = list()
datrain[[1]] = list(vname="precl",desc="PREC-L"       ,colour="firebrick"  ,order=5)
datrain[[2]] = list(vname="udel" ,desc="UDel-3.01"    ,colour="darkorange" ,order=4)
datrain[[3]] = list(vname="gpcc" ,desc="GPCC-6.0"     ,colour="yellow3"    ,order=6)
datrain[[4]] = list(vname="sheff",desc="PGMF"         ,colour="steelblue"  ,order=3)
datrain[[5]] = list(vname="gpcp" ,desc="GPCP-2.2"     ,colour="deepskyblue",order=2)
datrain[[6]] = list(vname="trmm" ,desc="TRMM/3B43-7.0",colour="purple3"    ,order=1)
#------------------------------------------------------------------------------------------#




#----- ED-2.2 simulation variables. -------------------------------------------------------#
ed22var       = list()
ed22var[[ 1]] = list( vname   = "et" 
                    , desc    = "Evapotranspiration"
                    , unit    = "mmoyr"
                    , plog    = FALSE
                    , min     = 600
                    , max     = 1600
                    , cscheme = "ipanoply"
                    )#end list
ed22var[[ 2]] = list( vname   = "transp"
                    , desc    = "Transpiration"
                    , unit    = "mmoyr"
                    , plog    = FALSE
                    , min     = 400
                    , max     = 1000
                    , cscheme = "clife"
                    )#end list
ed22var[[ 3]] = list( vname   = "gpp"
                    , desc    = "Gross Primary Productivity"
                    , unit    = "kgcom2"
                    , plog    = FALSE
                    , min     = 1.5
                    , max     = 3.5
                    , cscheme = "clife"
                    )#end list
ed22var[[ 4]] = list( vname   = "npp"
                    , desc    = "Net Primary Productivity"
                    , unit    = "kgcom2"
                    , plog    = FALSE
                    , min     = 1.2
                    , max     = 3.0
                    , cscheme = "clife"
                    )#end list
ed22var[[ 5]] = list( vname   = "rshort"
                    , desc    = "Incoming SW Radiation"
                    , unit    = "wom2"
                    , plog    = FALSE
                    , min     = 170
                    , max     = 230
                    , cscheme = "icloudy"
                    )#end list
ed22var[[ 6]] = list( vname   = "atm.temp"
                    , desc    = "Air temperature"
                    , unit    = "degc"
                    , plog    = FALSE
                    , min     = 15
                    , max     = 27
                    , cscheme = "panoply"
                    )#end list
ed22var[[ 7]] = list( vname   = "atm.vpd"
                    , desc    = "Air VPD"
                    , unit    = "hpa"
                    , plog    = FALSE
                    , min     = 5
                    , max     = 21
                    , cscheme = "panoply"
                    )#end list
ed22var[[ 8]] = list( vname   = "agb"
                    , desc    = "Above-ground biomass"
                    , unit    = "kgcom2"
                    , plog    = FALSE
                    , min     = 0
                    , max     = 22
                    , cscheme = "clife"
                    )#end list
ed22var[[ 9]] = list( vname   = "bsa"
                    , desc    = "Basal Area"
                    , unit    = "cm2om2"
                    , plog    = FALSE
                    , min     = 0
                    , max     = 32
                    , cscheme = "clife"
                    )#end list
ed22var[[10]] = list( vname   = "lai"
                    , desc    = "Leaf area index"
                    , unit    = "m2lom2"
                    , plog    = FALSE
                    , min     = 1.5
                    , max     = 4.5
                    , cscheme = "clife"
                    )#end list
#------------------------------------------------------------------------------------------#



#----- Statistics. ------------------------------------------------------------------------#
statvar       = list()
statvar[[ 1]] = list( vname   = "location"
                    , desc    = "Location Parameter"
                    , cscheme = "ipanoply"
                    , unit    = "mmoyr"
                    , min     = 900
                    , max     = 3500
                    , plog    = FALSE
                    )#end list
statvar[[ 2]] = list( vname   = "scale"
                    , desc    = "Scale Parameter"
                    , cscheme = "hue.hot"
                    , unit    = "mmoyr"
                    , min     = 100
                    , max     = 1000
                    , plog    = FALSE
                    )#end list
statvar[[ 3]] = list( vname   = "shape"
                    , desc    = "Shape Parameter"
                    , cscheme = "ipanoply"
                    , unit    = "empty"
                    , min     = -5
                    , max     = +5
                    , plog    = FALSE
                    )#end list
statvar[[ 4]] = list( vname   = "dlocation.rel"
                    , desc    = "Vulnerability (Location)"
                    , cscheme = "hue.hot"
                    , unit    = "empty"
                    , min     = -3
                    , max     =  0
                    , plog    = FALSE
                    )#end list
statvar[[ 5]] = list( vname   = "dscale.rel"
                    , desc    = "Vulnerability (Scale)"
                    , cscheme = "ihue.hot"
                    , unit    = "empty"
                    , min     =  0.0
                    , max     =  3.0
                    , plog    = FALSE
                    )#end list
statvar[[ 6]] = list( vname   = "dshape.abs"
                    , desc    = "Vulnerability (Shape)"
                    , cscheme = "hue.hot"
                    , unit    = "empty"
                    , min     = -3.0
                    , max     =  0.0
                    , plog    = FALSE
                    )#end list
statvar[[ 7]] = list( vname   = "mean"
                    , desc    = "Mean"
                    , cscheme = "ipanoply"
                    , unit    = "mmoyr"
                    , min     = 900
                    , max     = 3500
                    , plog    = FALSE
                    )#end list
statvar[[ 8]] = list( vname   = "sdev"
                    , desc    = "Standard Deviation"
                    , cscheme = "panoply"
                    , unit    = "mmoyr"
                    , min     = 100
                    , max     = 800
                    , plog    = FALSE
                    )#end list
statvar[[ 9]] = list( vname   = "skew"
                    , desc    = "Skewness"
                    , cscheme = "ipanoply"
                    , unit    = "empty"
                    , min     = -1.5
                    , max     =  1.5
                    , plog    = FALSE
                    )#end list
statvar[[10]] = list( vname   = "sn.cvar"
                    , desc    = "Coeff. of Variation (SN)"
                    , cscheme = "panoply"
                    , unit    = "empty"
                    , min     = 0.1
                    , max     = 0.3
                    , plog    = FALSE
                    )#end list
statvar[[11]] = list( vname   = "norm.cvar"
                    , desc    = "Coeff. of Variation"
                    , cscheme = "panoply"
                    , unit    = "empty"
                    , min     = 0.10
                    , max     = 0.25
                    , plog    = FALSE
                    )#end list
#------------------------------------------------------------------------------------------#




#----- Variables to be the colours for the 3-D plot. --------------------------------------#
three.dim      = list() 
three.dim[[1]] = list( vname = "agb", key = "AGB", unit = "kgcom2")
three.dim[[2]] = list( vname = "lai", key = "LAI", unit = "m2lom2")
#------------------------------------------------------------------------------------------#



#----- Years to use for Sheffield. --------------------------------------------------------#
sheff.yeara  = 1969
sheff.yearz  = 2008
#------------------------------------------------------------------------------------------#


#----- Texture-dependent return period. ---------------------------------------------------#
stext.pret  = c(    2,   6,   8,  16, 11)
lookup.pret = c(1.740,1.95,2.26,4.73,7.6)
use.stext   = c(2,2,2,8,8,6,8,8,6,8,11,16,16,8,11,16,11)
#------------------------------------------------------------------------------------------#


#------ List of soil types. ---------------------------------------------------------------#
stext.leg = c(seq(from=1,to=11,by=1),seq(from=14,to=17,by=1))
#------------------------------------------------------------------------------------------#


#----- Climate change settings. -----------------------------------------------------------#
n.years = 1
nsample = 10000
#------------------------------------------------------------------------------------------#



#----- List of points to be plotted. ------------------------------------------------------#
site       = list()
site[[ 1]] = list( iata = "gyf"
                 , lon  = -52.91
                 , lat  =   5.28
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 2]] = list( iata = "s67"
                 , lon  = -54.96
                 , lat  =  -2.86
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 3]] = list( iata = "bvb"
                 , lon  = -60.61
                 , lat  =   2.92
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 4]] = list( iata = "sip"
                 , lon  = -55.90
                 , lat  =   1.90
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 5]] = list( iata = "chb"
                 , lon  = -55.00
                 , lat  =  -9.00
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 6]] = list( iata = "scz"
                 , lon  = -62.50
                 , lat  = -14.50
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 7]] = list( iata = "pnd"
                 , lon  = -66.50
                 , lat  = -10.50
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 8]] = list( iata = "pcl"
                 , lon  = -74.57
                 , lat  =  -8.38
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
site[[ 9]] = list( iata = "mab"
                 , lon  = -49.14
                 , lat  =  -5.37
                 , col1 = "deeppink3"
                 , col2 = "turquoise2"
                 , fnt  = 2
                 , cex  = 0.6
                 )#end list
#------------------------------------------------------------------------------------------#



#----- List of interesting points for future reference. -----------------------------------#
ptref       = list()
ptref[[ 1]] = list(lab="Roraima"    ,lon=-61.1,lat=  1.2,col="white"    ,fnt=1,cex=0.6 )
ptref[[ 2]] = list(lab="Amazonas"   ,lon=-63.5,lat= -3.9,col="white"    ,fnt=1,cex=0.7 )
ptref[[ 3]] = list(lab="Pará"       ,lon=-52.7,lat= -4.3,col="white"    ,fnt=1,cex=0.7 )
ptref[[ 4]] = list(lab="Amapá"      ,lon=-51.9,lat=  1.5,col="white"    ,fnt=1,cex=0.6 )
ptref[[ 5]] = list(lab="Maranhão"   ,lon=-46.0,lat= -4.9,col="white"    ,fnt=1,cex=0.55)
ptref[[ 6]] = list(lab="Tocantins"  ,lon=-48.5,lat=-10.3,col="white"    ,fnt=1,cex=0.55)
ptref[[ 7]] = list(lab="Mato Grosso",lon=-55.4,lat=-12.8,col="white"    ,fnt=1,cex=0.7 )
ptref[[ 8]] = list(lab="Rondônia"   ,lon=-62.6,lat=-11.5,col="white"    ,fnt=1,cex=0.6 )
ptref[[ 9]] = list(lab="Acre"       ,lon=-71.0,lat= -8.7,col="white"    ,fnt=1,cex=0.7 )
ptref[[10]] = list(lab="Brazil"     ,lon=-58.1,lat= -8.0,col="black"    ,fnt=2,cex=0.9 )
ptref[[11]] = list(lab="Bolivia"    ,lon=-64.7,lat=-15.0,col="black"    ,fnt=2,cex=0.8 )
ptref[[12]] = list(lab="Peru"       ,lon=-76.5,lat= -6.3,col="black"    ,fnt=2,cex=0.9 )
ptref[[13]] = list(lab="Ecuador"    ,lon=-78.5,lat= -1.0,col="black"    ,fnt=2,cex=0.6 )
ptref[[14]] = list(lab="Colombia"   ,lon=-72.7,lat=  1.0,col="black"    ,fnt=2,cex=0.7 )
ptref[[15]] = list(lab="Venezuela"  ,lon=-64.5,lat=  7.0,col="black"    ,fnt=2,cex=0.7 )
ptref[[16]] = list(lab="GUY"        ,lon=-59.0,lat=  5.2,col="black"    ,fnt=2,cex=0.6 )
ptref[[17]] = list(lab="SUR"        ,lon=-56.1,lat=  4.2,col="black"    ,fnt=2,cex=0.6 )
ptref[[18]] = list(lab="GUF"        ,lon=-53.2,lat=  3.5,col="black"    ,fnt=2,cex=0.6 )
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#----- Load libraries and functions. ------------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#


#----- List of variables. -----------------------------------------------------------------#
nsites = length(site)
nptref = length(ptref)
#------------------------------------------------------------------------------------------#


#------ List of schemes to use col1. ------------------------------------------------------#
col1.schemes = c("panoply","ipanoply","muitas","imuitas","hue.purple","ihue.purple"
                ,"cloudy","icloudy","hue.blue","ihue.blue","hue.green","ihue.green"
                ,"clife","iclife","hue.cold","ihue.cold")
#------------------------------------------------------------------------------------------#




#----- Get dimensions. --------------------------------------------------------------------#
ndatrain   = length(datrain)
ned22var   = length(ed22var)
nstatvar   = length(statvar)
nthree.dim = length(three.dim)
#------------------------------------------------------------------------------------------#


#----- Update the reference table for return period based on soil types. ------------------#
idx.stext      = match(use.stext,stext.pret)
ref.table.pret = lookup.pret[idx.stext]
#------------------------------------------------------------------------------------------#


#----- Number of testable soil textures. --------------------------------------------------#
n.stext.pret = length(stext.pret)
#------------------------------------------------------------------------------------------#



#----- Standardise output format. ---------------------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#


#----- Create output paths in case they don't exist. --------------------------------------#
if (! file.exists(outroot   )) dir.create(outroot   )
if (! file.exists(rdata.path)) dir.create(rdata.path)
#------------------------------------------------------------------------------------------#



#----- Find out how many columns and rows to be plotted in multiple-panel plots. ----------#
lo.datrain = pretty.box(ndatrain)
#------------------------------------------------------------------------------------------#



#----- Find properties of the default soil types. -----------------------------------------#
cat(" + Find properties of the default soil types...","\n")
n.soil.ed21     = length(stext.leg)
soil.ed21       = sapply( X   = data.frame( t( mapply( FUN      = soil.params
                                                     , ntext    = stext.leg
                                                     , MoreArgs = list(isoilflg =  1
                                                                      ,slxsand  = -1
                                                                      ,slxclay  = -1
                                                                      )#end list
                                                     )#end mapply
                                            )#end t
                                          , stringsAsFactors = FALSE
                                        )#end data.frame
                        , FUN = unlist
                        )#end sapply
soil.ed21        = data.frame(soil.ed21,stringsAsFactors=FALSE)
for (n in names(soil.ed21)[! names(soil.ed21) %in% c("name","key")]){
   soil.ed21[[n]] = as.numeric(soil.ed21[[n]])
}#end for
soil.ed21$colour = stext.cols[soil.ed21$ntext]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Transform datrain into a data frame.                                               #
#------------------------------------------------------------------------------------------#
datrain = data.frame( apply( X = sapply( X = datrain, FUN = c), MARGIN = 1, FUN = unlist)
                    , stringsAsFactors = FALSE
                    )#end data.frame
for (iu in c("order")) datrain[[iu]] = as.numeric(datrain[[iu]])
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Transform site into a data frame.                                                  #
#------------------------------------------------------------------------------------------#
site   = data.frame( apply( X = sapply( X = site, FUN = c), MARGIN = 1, FUN = unlist)
                   , stringsAsFactors = FALSE
                   )#end data.frame
for (iu in c("lon","lat","fnt","cex")) site[[iu]] = as.numeric(site[[iu]])
#----- Fix the colours in case of a dark background. --------------------------------------#
if (ibackground %in% c(1,2)){
   black = site$col1 == "black"
   white = site$col1 == "white"
   site$col1[black]   = "white"
   site$col1[white]   = "black"
   black = site$col2 == "black"
   white = site$col2 == "white"
   site$col2[black]   = "white"
   site$col2[white]   = "black"
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Transform reference labels into a data frame.                                      #
#------------------------------------------------------------------------------------------#
ptref = data.frame( apply( X = sapply( X = ptref, FUN = c), MARGIN = 1, FUN = unlist)
                  , stringsAsFactors = FALSE
                  )#end data.frame
for (iu in c("lon","lat","cex","fnt")) ptref[[iu]] = as.numeric(ptref[[iu]])
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether to calculate or read in vulnerability maps.                            #
#------------------------------------------------------------------------------------------#
if (retrieve.data & file.exists(rdata.vulnerable)){
   #----- Load previous data. -------------------------------------------------------------#
   cat(" + Reloading vulnerability data from ",basename(rdata.vulnerable),"...","\n")
   load(rdata.vulnerable)
   nregion = nrow(region)
   nyear   = length(ymean$year)
   #---------------------------------------------------------------------------------------#


   #----- Find limits for plot. -----------------------------------------------------------#
   limlon  = range(region$lon)
   limlat  = range(region$lat)
   wlimlon = limlon + c(-1,1)  * diff(limlon) * lo.datrain$ncol
   wlimlat = limlat + c(-1,1)  * diff(limlat) * lo.datrain$nrow
   slimlon = limlon + c(0,0.3) * diff(limlon)
   slimlat = limlat
   size   = plotsize(proje=TRUE ,limlon=limlon ,limlat=limlat ,extendfc=TRUE ,paper=paper)
   ssize  = plotsize(proje=TRUE ,limlon=slimlon,limlat=slimlat,extendfc=FALSE,paper=paper)
   wsize  = plotsize(proje=TRUE ,limlon=wlimlon,limlat=wlimlat,extendfc=TRUE ,paper=paper)
   fsize  = plotsize(proje=FALSE,paper=paper)
   #---------------------------------------------------------------------------------------#

}else{

   #----- Find data again. ----------------------------------------------------------------#
   cat(" + Calculating vulnerability...","\n")
   #---------------------------------------------------------------------------------------#


   #----- Read in the statistics for the gridded data. ------------------------------------#
   cat ("   - Loading gridded data...","\n")
   load(rdata.gridded)
   #---------------------------------------------------------------------------------------#




   #----- First we find which stations we will keep. --------------------------------------#
   wmo$nyears = apply( X = is.finite(wmo$obs), MARGIN = 1, FUN = sum )
   use        = wmo$nyears >= wmo.nyears.min
   wmo$lon    = wmo$lon   [use  ]
   wmo$lat    = wmo$lat   [use  ]
   wmo$lola   = wmo$lola  [use, ]
   wmo$obs    = wmo$obs   [use, ]
   wmo$mod    = wmo$mod   [use,,]
   wmo$res    = wmo$res   [use,,]
   wmo$nyears = wmo$nyears[use  ]
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Find the site statistics.                                                        #
   #---------------------------------------------------------------------------------------#
   wmo.stats       = apply(X=wmo$obs,MARGIN=1,FUN=sn.stats,na.rm=TRUE)
   wmo$location    = wmo.stats["location",]
   wmo$scale       = wmo.stats["scale"   ,]
   wmo$shape       = wmo.stats["shape"   ,]
   wmo$mean        = apply(X=wmo$obs,MARGIN=1,FUN=mean    ,na.rm=TRUE)
   wmo$sdev        = apply(X=wmo$obs,MARGIN=1,FUN=sd      ,na.rm=TRUE)
   wmo$skew        = apply(X=wmo$obs,MARGIN=1,FUN=skew    ,na.rm=TRUE)
   wmo$location    = array(data=wmo$location,dim=dim(wmo$obs),dimnames=dimnames(wmo$obs))
   wmo$scale       = array(data=wmo$scale   ,dim=dim(wmo$obs),dimnames=dimnames(wmo$obs))
   wmo$shape       = array(data=wmo$shape   ,dim=dim(wmo$obs),dimnames=dimnames(wmo$obs))
   wmo$mean        = array(data=wmo$mean    ,dim=dim(wmo$obs),dimnames=dimnames(wmo$obs))
   wmo$sdev        = array(data=wmo$sdev    ,dim=dim(wmo$obs),dimnames=dimnames(wmo$obs))
   wmo$skew        = array(data=wmo$skew    ,dim=dim(wmo$obs),dimnames=dimnames(wmo$obs))
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      Remove information from data sets we are not going to use here.                  #
   #---------------------------------------------------------------------------------------#
   empty   = NA + wmo$obs
   dnames  = modifyList( x = dimnames(wmo$mod)
                       , val = list(model=c(dimnames(wmo$mod)$model,"sheff"))
                       )#end modifyList
   wmo$mod = abind(wmo$mod,empty,along=3,new.names=dnames)
   wmo$res = abind(wmo$res,empty,along=3,new.names=dnames)
   #---------------------------------------------------------------------------------------#






   #----- Initialise some variables to be read from the model output. ---------------------#
   nregion = nrow(region)
   region$stext          = rep(NA,times=nregion)
   region$et             = rep(NA,times=nregion)
   region$transp         = rep(NA,times=nregion)
   region$gpp            = rep(NA,times=nregion)
   region$npp            = rep(NA,times=nregion)
   region$rshort         = rep(NA,times=nregion)
   region$atm.temp       = rep(NA,times=nregion)
   region$atm.vpd        = rep(NA,times=nregion)
   region$agb            = rep(NA,times=nregion)
   region$bsa            = rep(NA,times=nregion)
   region$lai            = rep(NA,times=nregion)
   region$sheff.location = rep(NA,times=nregion)
   region$sheff.scale    = rep(NA,times=nregion)
   region$sheff.shape    = rep(NA,times=nregion)
   region$sheff.mean     = rep(NA,times=nregion)
   region$sheff.sdev     = rep(NA,times=nregion)
   region$sheff.skew     = rep(NA,times=nregion)
   region$sheff.bias     = rep(NA,times=nregion)
   region$sheff.sigres   = rep(NA,times=nregion)
   region$sheff.skew     = rep(NA,times=nregion)
   region$sheff.yeara    = rep(NA,times=nregion)
   region$sheff.yearz    = rep(NA,times=nregion)
   #---------------------------------------------------------------------------------------#


   #----- Find limits for plot. -----------------------------------------------------------#
   limlon  = range(region$lon)
   limlat  = range(region$lat)
   wlimlon = limlon + c(-1,1)  * diff(limlon) * lo.datrain$ncol
   wlimlat = limlat + c(-1,1)  * diff(limlat) * lo.datrain$nrow
   slimlon = limlon + c(0,0.3) * diff(limlon)
   slimlat = limlat
   size    = plotsize(proje=TRUE ,limlon=limlon ,limlat=limlat ,extendfc=TRUE ,paper=paper)
   ssize   = plotsize(proje=TRUE ,limlon=slimlon,limlat=slimlat,extendfc=FALSE,paper=paper)
   wsize   = plotsize(proje=TRUE ,limlon=wlimlon,limlat=wlimlat,extendfc=TRUE ,paper=paper)
   fsize   = plotsize(proje=FALSE,paper=paper)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over the rainfall methods, and define two variables: one for the absolute    #
   # change in local rainfall and another for the relative change in rainfall.             #
   #---------------------------------------------------------------------------------------#
   for (d in sequence(ndatrain)){
      region[[paste(datrain$vname[d],"dlocation","abs" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dlocation","rel" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dscale"   ,"abs" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dscale"   ,"rel" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dshape"   ,"abs" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dmean"    ,"abs" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dmean"    ,"rel" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dsdev"    ,"abs" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dsdev"    ,"rel" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"dskew"    ,"abs" ,sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"norm"     ,"cvar",sep=".")]] = rep(NA,times=nregion)
      region[[paste(datrain$vname[d],"sn"       ,"cvar",sep=".")]] = rep(NA,times=nregion)
   }#end for (d in sequence(ndatrain))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Save Sheffield rainfall to estimate the error.                                   #
   #---------------------------------------------------------------------------------------#
   sheff.year  = seq(from=sheff.yeara,to=sheff.yearz,by=1)
   sheff.yrain = matrix( data     = NA
                       , nrow     = nregion
                       , ncol     = length(sheff.year)
                       , dimnames = list(region$iata, sheff.year)
                       )#end matrix
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Loop over every polygon and find the ED-2.1 run.                                 #
   #---------------------------------------------------------------------------------------#
   for (r in sequence(nregion)){
      polyname    = region$pname[r]
      polyiata    = region$iata [r]
      polylon     = region$lon  [r]
      polylat     = region$lat  [r]


      #----- File info for this polygon. --------------------------------------------------#
      ed22.path   = file.path(here,polyname,"rdata_simple")
      ed22.status = file.path(ed22.path,paste("status_",polyname,".txt",sep=""))
      ed22.rdata  = file.path(ed22.path,paste(polyname,".RData",sep=""))
      #------------------------------------------------------------------------------------#

      #----- We will include only those polygons that finished. ---------------------------#
      cestfini    = file.exists(ed22.rdata) && file.exists(ed22.status)
      if (cestfini){
         cestfini = all(c(read.table(file=ed22.status)) == rdata.fini)
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Calculate vulnerability if data are available.                                 #
      #------------------------------------------------------------------------------------#
      if (cestfini){
         #----- Entertain the user. -------------------------------------------------------#
         cat ("     * Processing polygon: ",r,"/",nregion," (",polyiata,")...","\n")
         #---------------------------------------------------------------------------------#



         #----- Load polygon data. --------------------------------------------------------#
         load(ed22.rdata)
         emean = datum$emean
         #---------------------------------------------------------------------------------#


         #----- Select the period with useful data. ---------------------------------------#
         sel                   = datum$year >= sheff.yeara & datum$year <= sheff.yearz
         when                  = datum$when [sel]
         month                 = datum$month[sel]
         year                  = datum$year [sel]
         yeara.now             = min(year)
         yearz.now             = max(year)
         region$sheff.yeara[r] = yeara.now
         region$sheff.yearz[r] = yearz.now
         #---------------------------------------------------------------------------------#


         #----- Initiliase list with time series in case it isn't here. -------------------#
         if (! "ymean" %in% ls()){
            un.year = sort(unique(year))
            nyear   = length(unique(year))
            empty   = array( data     = NA
                           , dim      = c(nyear,nregion)
                           , dimnames = list(un.year,region$iata)
                           )#end array
            ymean   = list( year     = year
                          , lai      = empty
                          , agb      = empty
                          , bsa      = empty
                          , dlen     = empty
                          , et       = empty
                          , transp   = empty
                          , rain     = empty
                          , atm.temp = empty
                          , atm.vpd  = empty
                          , rshort   = empty
                          , stext    = empty
                          )#end list
            drought = list()
            for (s in sequence(n.stext.pret)){
               stext.key    = paste("stext",sprintf("%2.2i",stext.pret[s]),sep="")
               drought[[stext.key]] = list( lai      = NULL
                                          , agb      = NULL
                                          , bsa      = NULL
                                          , dlen     = NULL
                                          , et       = NULL
                                          , transp   = NULL
                                          , rain     = NULL
                                          , atm.temp = NULL
                                          , atm.vpd  = NULL
                                          , rshort   = NULL
                                          , et0      = NULL
                                          )#end list
            }#end for
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Grab a few average value variables from the data set. ---------------------#
         region$et        [r] = mean(emean$wflxca  [sel] * yr.day)
         region$transp    [r] = mean(emean$transp  [sel] * yr.day)
         region$gpp       [r] = mean(emean$gpp     [sel])
         region$npp       [r] = mean(emean$npp     [sel])
         region$rshort    [r] = mean(emean$rshort  [sel])
         region$atm.temp  [r] = mean(emean$atm.temp[sel])
         region$atm.vpd   [r] = mean(emean$atm.vpd [sel])
         region$agb       [r] = mean(emean$agb     [sel])
         region$lai       [r] = mean(emean$lai     [sel])
         region$bsa       [r] = mean(emean$basarea [sel])
         region$stext     [r] = datum$ntext
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Recalculate the water deficit.                                               #
         #---------------------------------------------------------------------------------#
         water.deficit    = rep(NA,times=datum$ntimes)
         nmon.wdef        = rep(NA,times=datum$ntimes)
         count            = rep(NA,times=datum$ntimes)
         mondays          = daymax(datum$when)
         water.deficit[1] = max(0, emean$wflxca[1]*mondays[1] - emean$rain[1])
         nmon.wdef    [1] = as.numeric(water.deficit[1] > 10)
         count        [1] = ifelse(nmon.wdef[1] == 0,NA,1)
         for (m in sequence(datum$ntimes)[-1]){
            water.deficit[m] = max(0, water.deficit[m-1]
                                    + datum$emean$wflxca[m]*mondays[m] - emean$rain[m])
            nmon.wdef    [m] = as.numeric(water.deficit[m] > 10) * (nmon.wdef[m-1] + 1)
            count        [m] = ifelse( nmon.wdef[m] == 0
                                     , NA
                                     , ifelse( nmon.wdef[m] == 1
                                             , max(c(0,count),na.rm=TRUE)+1
                                             , count[m-1]
                                             )#end ifelse
                                     )#end ifelse
         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Select the 12-month period by the end of the drought, or the entire        #
         # drought length in case it is a long drought.                                    #
         #---------------------------------------------------------------------------------#
         dry       = (! is.na(count)) & sel
         if (any(dry)){
            dry.dlen  = tapply(X=nmon.wdef[dry],INDEX=count[dry],FUN=max)


            idx       = sequence(datum$ntimes)
            last      = tapply(X=idx[dry],INDEX=count[dry],FUN=max)
            first     = pmax(1,pmin(tapply(X=idx[dry],INDEX=count[dry],FUN=min),last-11))
            dry.span  = floor( unlist( mapply( FUN      = seq
                                             , from     = first
                                             , to       = last
                                             , MoreArgs = list(by=1)
                                             ) ) )
            dry.count = last-first+1
            seq.dry   = seq_along(unique(count[dry]))
            dry.idx   = unlist(mapply(FUN=rep,x=seq.dry,each=dry.count))
            dry       = dry.span



            dry.lai      = tapply(emean$lai      [dry]             ,dry.idx,mean)
            dry.agb      = tapply(emean$agb      [dry]             ,dry.idx,mean)
            dry.bsa      = tapply(emean$basarea  [dry]             ,dry.idx,mean)
            dry.et       = tapply(emean$et       [dry]             ,dry.idx,mean)*yr.day
            dry.transp   = tapply(emean$transp   [dry]             ,dry.idx,mean)
            dry.rain     = tapply(emean$rain     [dry]/mondays[dry],dry.idx,mean)*yr.day
            dry.atm.temp = tapply(emean$atm.temp [dry]             ,dry.idx,mean)
            dry.atm.vpd  = tapply(emean$atm.vpd  [dry]             ,dry.idx,mean)
            dry.rshort   = tapply(emean$rshort   [dry]             ,dry.idx,mean)
            dry.et0      = region$et[r] + 0 * dry.lai
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Append the data to the lists. ---------------------------------------------#
         stext.key            = paste( "stext"
                                     , sprintf("%2.2i",use.stext[datum$ntext])
                                     , sep = ""
                                     )#end paste
         drought[[stext.key]]$lai      = c(drought[[stext.key]]$lai     , dry.lai     )
         drought[[stext.key]]$agb      = c(drought[[stext.key]]$agb     , dry.agb     )
         drought[[stext.key]]$bsa      = c(drought[[stext.key]]$bsa     , dry.bsa     )
         drought[[stext.key]]$dlen     = c(drought[[stext.key]]$dlen    , dry.dlen    )
         drought[[stext.key]]$et       = c(drought[[stext.key]]$et      , dry.et      )
         drought[[stext.key]]$et0      = c(drought[[stext.key]]$et0     , dry.et0     )
         drought[[stext.key]]$transp   = c(drought[[stext.key]]$transp  , dry.transp  )
         drought[[stext.key]]$rain     = c(drought[[stext.key]]$rain    , dry.rain    )
         drought[[stext.key]]$atm.temp = c(drought[[stext.key]]$atm.temp, dry.atm.temp)
         drought[[stext.key]]$atm.vpd  = c(drought[[stext.key]]$atm.vpd , dry.atm.vpd )
         drought[[stext.key]]$rshort   = c(drought[[stext.key]]$rshort  , dry.rshort  )
         #---------------------------------------------------------------------------------#

         #----- Save the annual means using the maximum water deficit of a given year. ----#
         peak   = tapply(nmon.wdef[sel],year,which.max) + seq(from=0,by=12,length.out=nyear)
         last   = pmax(1,peak-11)

         l12    = floor(mapply(FUN=seq,from=last,to=peak,MoreArgs=list(length.out=12)))
         yr.idx = rep(seq(from=yeara.now,to=yearz.now,by=1),each=12)

         ymean$dlen    [,r] = nmon.wdef[sel][peak]
         ymean$lai     [,r] = tapply(emean$lai     [sel]         ,year,mean)
         ymean$agb     [,r] = tapply(emean$agb     [sel]         ,year,mean)
         ymean$bsa     [,r] = tapply(emean$basarea [sel]         ,year,mean)
         ymean$et      [,r] = tapply(emean$et      [sel] * yr.day,year,mean)
         ymean$transp  [,r] = tapply(emean$transp  [sel] * yr.day,year,mean)
         ymean$rain    [,r] = tapply(emean$rain    [sel]         ,year,sum )
         ymean$atm.temp[,r] = tapply(emean$atm.temp[sel]         ,year,mean)
         ymean$atm.vpd [,r] = tapply(emean$atm.vpd [sel]         ,year,mean)
         ymean$rshort  [,r] = tapply(emean$rshort  [sel]         ,year,mean)
         ymean$stext   [,r] = datum$ntext
         #---------------------------------------------------------------------------------#



         #----- Calculate annual rainfall. ------------------------------------------------#
         emean.rain                = emean$rain[sel]
         sheff.yrain[r,]           = tapply(X=emean.rain,INDEX=year,FUN=sum)
         yr.stats                  = sn.stats(x=sheff.yrain[r,],na.rm=TRUE)
         region$sheff.location [r] = yr.stats[1]
         region$sheff.scale    [r] = yr.stats[2]
         region$sheff.shape    [r] = yr.stats[3]
         region$sheff.mean     [r] = mean(sheff.yrain[r,],na.rm=TRUE)
         region$sheff.sdev     [r] = sd  (sheff.yrain[r,],na.rm=TRUE)
         region$sheff.skew     [r] = skew(sheff.yrain[r,],na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #----- Pick the return period. ---------------------------------------------------#
         pret.now = ref.table.pret[region$stext[r]]
         #---------------------------------------------------------------------------------#


         #----- Loop over all rainfall data. ----------------------------------------------#
         for (d in sequence(ndatrain)){
            thisvname         = datrain$vname[d]
            thisdesc          = datrain$desc [d]
            cat ("       ~ Vulnerability based on: ",thisdesc,"...","\n")

            #----- Statistics for the current rainfall data set. --------------------------#
            now.location =   region[[paste(thisvname,"location",sep=".")]][r]
            now.scale    =   region[[paste(thisvname,"scale"   ,sep=".")]][r]
            now.shape    =   region[[paste(thisvname,"shape"   ,sep=".")]][r]
            now.mean     =   region[[paste(thisvname,"mean"    ,sep=".")]][r]
            now.sdev     =   region[[paste(thisvname,"sdev"    ,sep=".")]][r]
            now.skew     =   region[[paste(thisvname,"skew"    ,sep=".")]][r]
            now.nyears   = ( region[[paste(thisvname,"yearz"   ,sep=".")]][r]
                           - region[[paste(thisvname,"yeara"   ,sep=".")]][r] ) + 1
            #------------------------------------------------------------------------------#



            #----- Coefficients of variation. ---------------------------------------------#
            now.norm.cvar          = now.sdev / now.mean
            now.sn.cvar            = ( diff(qsn( p        = pnorm(c(-0.5,0.5))
                                               , location = now.location
                                               , scale    = now.scale
                                               , shape    = now.shape ) )
                                     / qsn( p        = 0.5
                                          , location = now.location
                                          , scale    = now.scale
                                          , shape    = now.shape ) )
            norm.cvar              = paste(thisvname,"norm","cvar",sep=".")
            sn.cvar                = paste(thisvname,"sn"  ,"cvar",sep=".")
            region[[norm.cvar]][r] = now.norm.cvar
            region[[sn.cvar  ]][r] = now.sn.cvar
            #------------------------------------------------------------------------------#




            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
            #    Find the probability based on return period.  If using number of years,   #
            # we may use the standard probability distribution functions, otherwise we     #
            # must use a simple simulator to estimate probability.                         #
            #------------------------------------------------------------------------------#
            now.prob     = n.years / pret.now
            now.drought  = n.years * region$et[r]

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #     Shift in location parameter.                                             #
            #------------------------------------------------------------------------------#
            cat ("         > Location shift...","\n")
               ntry         = 0
               again        = TRUE
               while (again & ntry < ntry.max){
                  ntry      = ntry + 1
                  delta     = 10. * pmax(abs(now.shape),1) * now.scale
                  now.lower = n.years * (now.location - ntry * delta)
                  now.upper = n.years * (now.location + ntry * delta)

                  now.lower = min(now.drought,now.lower)
                  now.upper = max(now.drought,now.upper)
                  #----- Use try to ensure it doesn't crash. ------------------------------#
                  ans       = try( qsn.mult( z        = now.drought
                                           , n        = n.years
                                           , p        = now.prob
                                           , location = now.location
                                           , scale    = now.scale
                                           , shape    = now.shape
                                           , shift    = "location"
                                           , lower    = now.lower
                                           , upper    = now.upper
                                           , nsample  = nsample
                                           )#end qsn.mult
                                 , silent = TRUE
                                 )#end try
                  again = "try-error" %in% is(ans)
                  if (again && loud){
                     cat ("           ~ Failed, try again (",ntry,"/",ntry.max,")...","\n")
                  }#end if
                  #------------------------------------------------------------------------#
               }#end while
               #---------------------------------------------------------------------------#



               #----- Copy answer to the dataset. -----------------------------------------#
               if (! "try-error" %in% is(ans)){
                  dlocation.abs              = paste(thisvname,"dlocation","abs",sep=".")
                  dlocation.rel              = paste(thisvname,"dlocation","rel",sep=".")
                  dmean.abs                  = paste(thisvname,"dmean"    ,"abs",sep=".")
                  dmean.rel                  = paste(thisvname,"dmean"    ,"rel",sep=".")
                  region[[dlocation.abs]][r] = ans["dlocation"]
                  region[[dmean.abs    ]][r] = ans["dmean"    ]
                  region[[dlocation.rel]][r] = ans["dlocation"] / now.scale
                  region[[dmean.rel    ]][r] = ans["dmean"    ] / now.sdev
               }else{
                  cat ("           ~ Failed finding location vulnerability...","\n")
               }#end if
               #---------------------------------------------------------------------------#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #     Shift in scale parameter.                                                #
            #------------------------------------------------------------------------------#
            cat ("         > Scale shift...","\n")
               ntry         = 0
               again        = TRUE
               while (again & ntry < ntry.max){
                  ntry      = ntry + 1
                  now.lower = now.scale / 100^ntry
                  now.upper = now.scale * 100^ntry
                  #----- Use try to ensure it doesn't crash. ------------------------------#
                  ans       = try( qsn.mult( z        = now.drought
                                           , n        = n.years
                                           , p        = now.prob
                                           , location = now.location
                                           , scale    = now.scale
                                           , shape    = now.shape
                                           , shift    = "scale"
                                           , lower    = now.lower
                                           , upper    = now.upper
                                           , nsample  = nsample
                                           )#end qsn.mult
                                 , silent = TRUE
                                 )#end try
                  again = "try-error" %in% is(ans)
                  if (again && loud){
                     cat ("           ~ Failed, try again (",ntry,"/",ntry.max,")...","\n")
                  }#end if
                  #------------------------------------------------------------------------#
               }#end while
               #---------------------------------------------------------------------------#



               #----- Copy answer to the dataset. -----------------------------------------#
               if (! "try-error" %in% is(ans)){
                  dscale.abs                 = paste(thisvname,"dscale"   ,"abs",sep=".")
                  dscale.rel                 = paste(thisvname,"dscale"   ,"rel",sep=".")
                  dsdev.abs                  = paste(thisvname,"dsdev"    ,"abs",sep=".")
                  dsdev.rel                  = paste(thisvname,"dsdev"    ,"rel",sep=".")
                  region[[dscale.abs   ]][r] = ans["dscale"   ]
                  region[[dsdev.abs    ]][r] = ans["dsdev"    ]
                  region[[dscale.rel   ]][r] = ans["dscale"   ] / now.scale
                  region[[dsdev.rel    ]][r] = ans["dsdev"    ] / now.sdev
               }else{
                  cat ("           ~ Failed finding scale vulnerability...","\n")
               }#end if
               #---------------------------------------------------------------------------#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #     Shift in shape parameter.                                                #
            #------------------------------------------------------------------------------#
            cat ("         > Shape shift...","\n")
               ntry         = 0
               again        = TRUE
               while (again & ntry < ntry.max){
                  ntry      = ntry + 1
                  now.lower = - abs(now.shape) * 10^ntry
                  now.upper = + abs(now.shape) * 10^ntry
                  #----- Use try to ensure it doesn't crash. ------------------------------#
                  ans       = try( qsn.mult( z        = now.drought
                                           , n        = n.years
                                           , p        = now.prob
                                           , location = now.location
                                           , scale    = now.scale
                                           , shape    = now.shape
                                           , shift    = "shape"
                                           , lower    = now.lower
                                           , upper    = now.upper
                                           , nsample  = nsample
                                           )#end qsn.mult
                                 , silent = TRUE
                                 )#end try
                  again = "try-error" %in% is(ans)
                  if (again && loud){
                     cat ("           ~ Failed, try again (",ntry,"/",ntry.max,")...","\n")
                  }#end if
                  #------------------------------------------------------------------------#
               }#end while
               #---------------------------------------------------------------------------#



               #----- Copy answer to the dataset. -----------------------------------------#
               if (! "try-error" %in% is(ans)){
                  dshape.abs                 = paste(thisvname,"dshape"   ,"abs",sep=".")
                  dskew.abs                  = paste(thisvname,"dskew"    ,"abs",sep=".")
                  region[[dshape.abs   ]][r] = ans["dshape"   ]
                  region[[dskew.abs    ]][r] = ans["dskew"    ]
               }else{
                  cat ("           ~ Failed finding shape vulnerability...","\n")
               }#end if
               #---------------------------------------------------------------------------#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         }#end for
         #---------------------------------------------------------------------------------#


         #----- Clear data. ---------------------------------------------------------------#
         rm(datum,emean,dlocation.abs,dlocation.rel,dscale.abs,dscale.rel,dshape.abs
           ,dmean.abs,dmean.rel,dsdev.abs,dsdev.rel,dskew.abs,now.prob,now.drought
           ,now.lower,now.upper)
         #---------------------------------------------------------------------------------#

      }else{
         #----- Entertain the user. -------------------------------------------------------#
         cat ("     * Skipping polygon: ",r,"/",nregion," (not finished)...","\n")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Find the nearest neighbour for all stations.                                      #
   #---------------------------------------------------------------------------------------#
   region.lola               = as.points(data.frame(x=region$lon,y=region$lat))
   idx.lola                  = nearest.neighbour(wmo$lola,region.lola)
   idx.year                  = match(wmo$year,unique(year))
   ok.year                   = which(! is.na(idx.year))
   idx.year                  = idx.year[ok.year]
   wmo$mod[,ok.year,"sheff"] = sheff.yrain[idx.lola,idx.year]
   wmo$res[       ,,"sheff"] = wmo$mod[,,"sheff"] - wmo$obs
   now.bias                  = mean(c(wmo$res[,,"sheff"]),na.rm=TRUE)
   now.sigres                = sd  (c(wmo$res[,,"sheff"]),na.rm=TRUE)
   now.n                     = sum(! is.na(wmo$res[,,"sheff"]))
   now.rmse                  = sqrt(now.bias^2+(now.n-1)*now.sigres^2/now.n)
   region$sheff.bias         = rep(now.bias  ,times=nregion)
   region$sheff.sigres       = rep(now.sigres,times=nregion)
   region$sheff.rmse         = rep(now.rmse  ,times=nregion)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Find the support (log-likelihood) of all rainfall data sets using the             #
   # distribution of the residuals.                                                        #
   #---------------------------------------------------------------------------------------#
   wmo$mod.bias   = apply(X=wmo$res           ,MARGIN=c(1,3),FUN=mean,na.rm=TRUE)
   wmo$mod.sigma  = apply(X=wmo$res           ,MARGIN=c(1,3),FUN=sd  ,na.rm=TRUE)
   wmo$mod.n      = apply(X=is.finite(wmo$res),MARGIN=c(1,3),FUN=sum ,na.rm=TRUE)
   wmo$mod.mse    = wmo$mod.bias^2+(wmo$mod.n-1)*wmo$mod.sigma^2 / wmo$mod.n
   wmo$mod.weight = wmo$mod.n / wmo$mod.mse
   #----- Convert sparse data into NA. ----------------------------------------------------#
   sparse                = wmo$mod.n < 2
   wmo$mod.bias  [sparse] = NA
   wmo$mod.sigma [sparse] = NA
   wmo$mod.n     [sparse] = 0
   wmo$mod.mse   [sparse] = NA
   wmo$mod.weight[sparse] = 0
   #----- Convert everything to data frames. ----------------------------------------------#
   wmo$mod.bias   = as.data.frame(wmo$mod.bias  )
   wmo$mod.sigma  = as.data.frame(wmo$mod.sigma )
   wmo$mod.n      = as.data.frame(wmo$mod.n     )
   wmo$mod.mse    = as.data.frame(wmo$mod.mse   )
   wmo$mod.weight = as.data.frame(wmo$mod.weight)
   #----- Find the weight for each data set. ----------------------------------------------#
   wmo$mod.weight = mapply( FUN      = weighted.mean
                          , x        = wmo$mod.weight
                          , w        = wmo$mod.n
                          , MoreArgs = list(na.rm=TRUE)
                          )#end mapply
   wmo$mod.weight = wmo$mod.weight / sum(wmo$mod.weight)
   #---------------------------------------------------------------------------------------#



   #----- Load previous data. -------------------------------------------------------------#
   cat(" + Saving vulnerability data to ",basename(rdata.vulnerable),"...","\n")
   save(region,ymean,drought,wmo,file=rdata.vulnerable)
   #---------------------------------------------------------------------------------------#
}#end if
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Grab the variables that we will plot and make lists.                                #
#------------------------------------------------------------------------------------------#



#------ Limits for the longitude and latitude. --------------------------------------------#
longitude.limit = range(region$lon,na.rm=TRUE)
latitude.limit  = range(region$lat,na.rm=TRUE)
lonplot         = pretty.lonlat(x=longitude.limit,n=4,type="lon")
latplot         = pretty.lonlat(x=latitude.limit ,n=4,type="lat")
yeara           = median(region$sheff.yeara,na.rm=TRUE)
yearz           = median(region$sheff.yearz,na.rm=TRUE)
#------------------------------------------------------------------------------------------#

#------ List for plotting the maps. -------------------------------------------------------#
list.lon           = list()
list.lat           = list()
list.dlon          = list()
list.dlat          = list()
list.location      = list()
list.scale         = list()
list.shape         = list()
list.mean          = list()
list.sdev          = list()
list.skew          = list()
list.norm.cvar     = list()
list.sn.cvar       = list()
list.dlocation.rel = list()
list.dscale.rel    = list()
list.dshape.abs    = list()
list.dmean.rel     = list()
list.dsdev.rel     = list()
list.dskew.abs     = list()
list.sub           = list()
list.plot.after    = list()
list.x.axis        = list()
list.y.axis        = list()
for (d in sequence(ndatrain)){
   #----- Grab information about this rainfall data set. ----------------------------------#
   thisvname               = datrain$vname[d]
   thisdesc                = datrain$desc [d]
   list.lon          [[d]] = region$lon
   list.lat          [[d]] = region$lat
   list.dlon         [[d]] = median(diff(sort(unique(region$lon))),na.rm=TRUE)
   list.dlat         [[d]] = median(diff(sort(unique(region$lat))),na.rm=TRUE)
   list.location     [[d]] = region[[paste(thisvname,"location"        ,sep=".")]]
   list.scale        [[d]] = region[[paste(thisvname,"scale"           ,sep=".")]]
   list.shape        [[d]] = region[[paste(thisvname,"shape"           ,sep=".")]]
   list.mean         [[d]] = region[[paste(thisvname,"mean"            ,sep=".")]]
   list.sdev         [[d]] = region[[paste(thisvname,"sdev"            ,sep=".")]]
   list.skew         [[d]] = region[[paste(thisvname,"skew"            ,sep=".")]]
   list.norm.cvar    [[d]] = region[[paste(thisvname,"norm"     ,"cvar",sep=".")]]
   list.sn.cvar      [[d]] = region[[paste(thisvname,"sn"       ,"cvar",sep=".")]]
   list.dlocation.rel[[d]] = region[[paste(thisvname,"dlocation","rel" ,sep=".")]]
   list.dscale.rel   [[d]] = region[[paste(thisvname,"dscale"   ,"rel" ,sep=".")]]
   list.dshape.abs   [[d]] = region[[paste(thisvname,"dshape"   ,"abs" ,sep=".")]]
   list.dmean.rel    [[d]] = region[[paste(thisvname,"dmean"    ,"rel" ,sep=".")]]
   list.dsdev.rel    [[d]] = region[[paste(thisvname,"dsdev"    ,"rel" ,sep=".")]]
   list.dskew.abs    [[d]] = region[[paste(thisvname,"dskew"    ,"abs" ,sep=".")]]
   list.x.axis       [[d]] = list(side=1,at=lonplot$at,labels=lonplot$labels)
   list.y.axis       [[d]] = list(side=2,at=latplot$at,labels=latplot$labels,las=1)
   #---------------------------------------------------------------------------------------#

   #----- Plotting options. ---------------------------------------------------------------#
   yeara                = median(region[[paste(thisvname,"yeara",sep=".")]],na.rm=TRUE)
   yearz                = median(region[[paste(thisvname,"yearz",sep=".")]],na.rm=TRUE)
   list.sub       [[d]] = list( main     = paste(thisdesc," (",yeara,"-",yearz,")",sep="")
                              , cex.main = 0.8 * cex.main
                              , line     = 0.5
                              )#end list
   list.plot.after[[d]] = list( southammap = list(lwd=1,col=grey.fg   )
                              , amazonmap  = list(lwd=2,col=foreground)
                              )#end list
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Plot some maps with the properties that do not depend on the rainfall statistics.    #
#------------------------------------------------------------------------------------------#
cat (" + Plot the Amazon map of a few regional properties...","\n")
for ( v in sequence(ned22var)){
   thised22      = ed22var[[v]]
   this.vname    = thised22$vname
   this.desc     = thised22$desc
   this.key.main = desc.unit(desc=NULL,unit=untab[[thised22$unit]])
   this.cscheme  = get(thised22$cscheme)
   this.plog     = thised22$plog
   this.min      = thised22$min
   this.max      = thised22$max


   #---- Select colour for sites. ---------------------------------------------------------#
   if (thised22$cscheme %in% col1.schemes){
      col.key = "col1"
   }else{
      col.key = "col2"
   }#end if
   #---------------------------------------------------------------------------------------#

   cat ("   - ",this.desc,"...","\n")


   #----- Find the Z scale based on the preferences. --------------------------------------#
   this.region    = region[[this.vname]]
   if (! (is.null(this.min) || is.null(this.max))){
      this.region = pmax(this.min,pmin(this.max,this.region))
      this.limit  = c(this.min,this.max)
   }else{
      this.limit  = range(this.region,finite=TRUE)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Make sure all classes are included. ---------------------------------------------#
   if (this.plog){
      this.levels  = pretty.log(this.limit,n=n.colourbar)
      this.levels  = c( min(this.levels) / mean(ediff(this.levels))
                      , this.levels
                      , max(this.levels) * mean(ediff(this.levels))
                      )#end c
      this.nlevels = length(this.levels)
   }else{
      this.levels  = pretty(this.limit,n=n.colourbar)
      this.levels  = c( min(this.levels) - mean(diff(this.levels))
                      , this.levels
                      , max(this.levels) + mean(diff(this.levels))
                      )#end c
      this.nlevels = length(this.levels)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Loop over output format. --------------------------------------------------------#
   for (o in sequence(nout)){
      #------ Open file. ------------------------------------------------------------------#
      fichier = file.path(outroot,paste("amzmap_",this.vname,".",outform[o],sep=""))
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
      #------------------------------------------------------------------------------------#


      #----- Map. -------------------------------------------------------------------------#
      image.map( x              = region$lon
               , y              = region$lat
               , z              = this.region
               , colour.palette = this.cscheme
               , levels         = this.levels
               , nlevels        = this.nlevels
               , na.col         = miss.colour
               , xlim           = longitude.limit
               , ylim           = latitude.limit
               , zlim           = this.limit
               , x.axis.options = list(side=1,at=lonplot$at,labels=lonplot$labels)
               , y.axis.options = list(side=2,at=latplot$at,labels=latplot$labels,las=1)
               , main.title     = list( text = paste("Mean ",this.desc
                                                    ," (",yeara,"-",yearz,")",sep="")
                                      , font = 2
                                      , cex  = cex.main
                                      )#end list
               , key.title      = title( main     = this.key.main
                                       , cex.main = 0.8 * cex.main
                                       , line     = 1
                                       )#end title
               , key.log        = this.plog
               , plot.after     = list( abline     = list( v   = lonplot$at
                                                         , h   = latplot$at
                                                         , lwd = 1
                                                         , lty = "dotted"
                                                         , col = grid.colour
                                                         )#end list
                                      , southammap = list(lwd=1,col=grey.fg   )
                                      , amazonmap  = list(lwd=2,col=foreground)
                                      , text       = list( x      = site$lon
                                                         , y      = site$lat
                                                         , labels = site$iata
                                                         , col    = site[[col.key]]
                                                         , cex    = 1.5*site$cex
                                                         , font   = site$fnt
                                                         )#end list
                                      )#end list
               , f.key          = 1/6
               )#end image.map
      #------------------------------------------------------------------------------------#


      #----- Close the device. ------------------------------------------------------------#
      if (outform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
      bye = clean.tmp()
      #------------------------------------------------------------------------------------#
   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(ned22var))
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Plot the statistics for all maps.                                                    #
#------------------------------------------------------------------------------------------#
cat (" + Plot the statistics...","\n")
for (v in sequence(nstatvar)){
   this.stat     = statvar[[v]]
   this.vname    = this.stat$vname
   
   this.desc     = this.stat$desc
   this.key.main = desc.unit(desc=NULL,unit=untab[[this.stat$unit]])
   this.cscheme  = get(this.stat$cscheme)
   this.min      = this.stat$min
   this.max      = this.stat$max
   this.plog     = this.stat$plog


   #---- Select colour for sites. ---------------------------------------------------------#
   if (this.stat$cscheme %in% col1.schemes){
      col.key = "col1"
   }else{
      col.key = "col2"
   }#end if
   #---------------------------------------------------------------------------------------#

   
   cat ("   - ",this.desc,"...","\n")


   #------ Limits for the z key scale. ----------------------------------------------------#
   if (! ( is.null(this.min) || is.null(this.max) )){
      this.list    = lapply( X   = lapply( X   = get(paste("list",this.vname,sep="."))
                                         , FUN = pmin
                                         , this.max
                                         )#end lapply
                           , FUN = pmax
                           , this.min
                           )#end lapply
      this.limit   = c(this.min,this.max)
   }else{
      this.list    = get(paste("list",this.vname,sep="."))
      this.limit   = range(unlist(this.list),na.rm=TRUE)
   }#end if
   if (this.plog){
      this.levels  = pretty.log(this.limit,n=n.colourbar)
      this.levels  = c( min(this.levels) / mean(ediff(this.levels))
                      , this.levels
                      , max(this.levels) * mean(ediff(this.levels))
                      )#end c
      this.nlevels = length(this.levels)
   }else{
      this.levels  = pretty(this.limit,n=n.colourbar)
      this.levels  = c( min(this.levels) - mean(diff(this.levels))
                      , this.levels
                      , max(this.levels) + mean(diff(this.levels))
                      )#end c
      this.nlevels = length(this.levels)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Loop over output format. --------------------------------------------------------#
   for (o in sequence(nout)){
      #------ Open file. ------------------------------------------------------------------#
      fichier = file.path(outroot,paste("statmap_",this.vname,".",outform[o],sep=""))
      if(outform[o] == "x11"){
         X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
      }else if(outform[o] == "png"){
         png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
            ,pointsize=ptsz,res=depth)
      }else if(outform[o] == "eps"){
         postscript(file=fichier,width=wsize$width,height=wsize$height
                   ,pointsize=ptsz,paper=wsize$paper)
      }else if(outform[o] == "pdf"){
         pdf(file=fichier,onefile=FALSE
            ,width=wsize$width,height=wsize$height,pointsize=ptsz,paper=wsize$paper)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Map. -------------------------------------------------------------------------#
      image.map( x              = list.lon
               , y              = list.lat
               , z              = this.list
               , dx             = list.dlon
               , dy             = list.dlat
               , colour.palette = this.cscheme
               , levels         = this.levels
               , nlevels        = this.nlevels
               , na.col         = miss.colour
               , xlim           = longitude.limit
               , ylim           = latitude.limit
               , x.axis.options = list.x.axis
               , y.axis.options = list.y.axis
               , sub.options    = list.sub
               , main.title     = list( text = paste("Annual Rainfall",this.desc,sep=" - ")
                                      , font = 2
                                      , cex  = cex.main
                                      )#end list
               , key.title      = title( main     = this.key.main
                                       , cex.main = 0.8 * cex.main
                                       , line     = 1
                                       )#end title
               , key.log        = this.plog
               , plot.after     = list( abline     = list( v   = lonplot$at
                                                         , h   = latplot$at
                                                         , lwd = 1
                                                         , lty = "dotted"
                                                         , col = grid.colour
                                                         )#end list
                                      , southammap = list(lwd=1,col=grey.fg   )
                                      , amazonmap  = list(lwd=2,col=foreground)
                                      , text       = list( x      = site$lon
                                                         , y      = site$lat
                                                         , labels = site$iata
                                                         , col    = site[[col.key]]
                                                         , cex    = 1.2*site$cex
                                                         , font   = site$fnt
                                                         )#end list
                                      )#end list
               , matrix.plot    = TRUE
               , edge.axes      = TRUE
               , f.key          = 1/9
               )#end image.map
      #------------------------------------------------------------------------------------#


      #----- Close the device. ------------------------------------------------------------#
      if (outform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
      bye = clean.tmp()
      #------------------------------------------------------------------------------------#
   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(nstatvar))
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Plot the statistics for all maps.                                                    #
#------------------------------------------------------------------------------------------#
cat (" + Plot the RMSE for each site and data set...","\n")

   #---- Select colour for sites. ---------------------------------------------------------#
   if ("panoply" %in% col1.schemes){
      col.key = "col1"
   }else{
      col.key = "col2"
   }#end if
   #---------------------------------------------------------------------------------------#


   #---- Get the counts and the root mean square error. -----------------------------------#
   rmse.use        = ( sqrt(wmo$mod.bias[,datrain$vname]^2+wmo$mod.sigma[,datrain$vname]^2)
                     / rep(wmo$sdev[,1],times=length(datrain$vname)) )
   count.use       = wmo$mod.n[,datrain$vname]
   skip            = count.use < n.lsq.min
   rmse.use [skip] = NA
   count.use[skip] = NA
   lon.use   = wmo$lon + 0 * rmse.use
   lat.use   = wmo$lat + 0 * rmse.use
   #---------------------------------------------------------------------------------------#



   #------ Limits for the z key scale. ----------------------------------------------------#
   this.limit   = range(unlist(rmse.use),na.rm=TRUE)
   this.levels  = pretty(this.limit,n=n.colourbar)
   this.levels  = c( min(this.levels) - mean(diff(this.levels))
                   , this.levels
                   , max(this.levels) + mean(diff(this.levels))
                   )#end c
   this.nlevels = length(this.levels)
   #---------------------------------------------------------------------------------------#


   #----- Get the number of available years and find the size scale. ----------------------#
   nyears.min  = min(unlist(count.use),na.rm=TRUE)
   nyears.max  = max(unlist(count.use),na.rm=TRUE)
   nyears.leg  = sort(unique(c(nyears.min,nyears.max,pretty(unlist(wmo$mod.n),n=4))))
   keep        = nyears.leg %in% seq(from=nyears.min,to=nyears.max,by=1)
   nyears.leg  = nyears.leg[keep]
   mod.n.cex   = 0.7 + 1.3 * ( count.use  - nyears.min ) / ( nyears.max - nyears.min )
   cex.leg     = 0.7 + 1.3 * ( nyears.leg - nyears.min ) / ( nyears.max - nyears.min )
   #---------------------------------------------------------------------------------------#


   #----- Create lists. -------------------------------------------------------------------#
   bycol     = unlist(col(rmse.use))
   lon.list  = split(x=unlist(lon.use  ),f=bycol); names(lon.list)  = datrain$desc
   lat.list  = split(x=unlist(lat.use  ),f=bycol); names(lat.list)  = datrain$desc
   rmse.list = split(x=unlist(rmse.use ),f=bycol); names(rmse.list) = datrain$desc
   cex.list  = split(x=unlist(mod.n.cex),f=bycol); names(cex.list)  = datrain$desc
   #---------------------------------------------------------------------------------------#



   #----- Pack the legend options to a list. ----------------------------------------------#
   xyz.legend  = list( x      = "center"
                     , inset  = 0.0
                     , legend = nyears.leg
                     , col    = foreground
                     , pch    = pch.wmo
                     , pt.cex = cex.leg
                     , title  = expression(bold("Comparable years"))
                     , ncol   = length(nyears.leg)
                     , cex    = 0.9 * cex.ptsz
                     , xpd    = TRUE
                     )#end legend
   #---------------------------------------------------------------------------------------#






   #----- Loop over output format. --------------------------------------------------------#
   for (o in sequence(nout)){
      #------ Open file. ------------------------------------------------------------------#
      fichier = file.path(outroot,paste("rmse_rain",outform[o],sep="."))
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
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Plot the skew plot.                                                           #
      #------------------------------------------------------------------------------------#
      par(par.user)
      xyz.plot  ( x              = lon.list
                , y              = lat.list
                , z              = rmse.list
                , fixed.xlim     = TRUE
                , fixed.ylim     = TRUE
                , xlim           = longitude.limit
                , ylim           = latitude.limit
                , colour.palette = panoply
                , levels         = this.levels
                , nlevels        = this.nlevels
                , pch            = pch.wmo
                , cex            = cex.list
                , xyz.before     = list( southammap = list(col=grey.mg)
                                       , amazonmap  = list(col=foreground,lwd=2.0)
                                       )#end list
                , xyz.after      = list( text       = list( x      = site$lon
                                                          , y      = site$lat
                                                          , labels = site$iata
                                                          , col    = site[[col.key]]
                                                          , cex    = 1.2*site$cex
                                                          , font   = site$fnt
                                                          )#end list
                                       )#end list
                , xyz.title      = list ( main = "Relative root mean square error")
                , key.title      = title( main = desc.unit(desc=NULL,unit=untab$empty))
                , xyz.xaxis      = list(side=1,at=lonplot$at,labels=lonplot$labels)
                , xyz.yaxis      = list(side=2,at=latplot$at,labels=latplot$labels,las=1)
                , xyz.legend     = xyz.legend
                , key.width      = 7.0
                , leg.height     = 8.0
                )#end skill.plot
      #------------------------------------------------------------------------------------#


      #----- Close the device. ------------------------------------------------------------#
      if (outform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
      bye = clean.tmp()
      #------------------------------------------------------------------------------------#
   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Plot the location for all maps.                                                      #
#------------------------------------------------------------------------------------------#
cat (" + Plot the statistics scaled by Error...","\n")
for (v in sequence(nstatvar)){
   this.stat     = statvar[[v]]
   this.vname    = this.stat$vname
   this.desc     = this.stat$desc
   this.key.main = desc.unit(desc=NULL,unit=untab[[this.stat$unit]])
   this.cscheme  = get(this.stat$cscheme)
   this.min      = this.stat$min
   this.max      = this.stat$max
   this.plog     = this.stat$plog

   cat ("   - ",this.desc,"...","\n")



   #---- Select colour for sites. ---------------------------------------------------------#
   if (this.stat$cscheme %in% col1.schemes){
      col.key = "col1"
   }else{
      col.key = "col2"
   }#end if
   #---------------------------------------------------------------------------------------#





   #------ Find the weighted average of the property. -------------------------------------#
   stat.list           = region[,paste(datrain$vname,this.stat$vname,sep=".")]
   weight              = matrix( data     = wmo$mod.weight[datrain$vname]
                               , nrow     = nregion
                               , ncol     = ndatrain
                               , dimnames = dimnames(stat.list)
                               , byrow    = TRUE
                               )#end matrix
   weighted.stat       = ( rowSums(stat.list*weight      ,na.rm=TRUE) 
                         / rowSums(0 * stat.list + weight,na.rm=TRUE) )
   miss                = ! is.finite(weighted.stat)
   weighted.stat[miss] = NA
   #---------------------------------------------------------------------------------------#




   #------ Limits for the z key scale. ----------------------------------------------------#
   if (! ( is.null(this.min) || is.null(this.max) )){
      weighted.stat = pmax(this.min,pmin(this.max,weighted.stat))
      this.limit    = c(this.min,this.max)
   }else{
      this.limit   = range(weighted.stat,na.rm=TRUE)
   }#end if
   if (this.plog){
      this.levels  = pretty.log(this.limit,n=n.colourbar)
      this.levels  = c( min(this.levels) / mean(ediff(this.levels))
                      , this.levels
                      , max(this.levels) * mean(ediff(this.levels))
                      )#end c
      this.nlevels = length(this.levels)
   }else{
      this.levels  = pretty(this.limit,n=n.colourbar)
      this.levels  = c( min(this.levels) - mean(diff(this.levels))
                      , this.levels
                      , max(this.levels) + mean(diff(this.levels))
                      )#end c
      this.nlevels = length(this.levels)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Loop over output format. --------------------------------------------------------#
   for (o in sequence(nout)){
      #------ Open file. ------------------------------------------------------------------#
      fichier = file.path(outroot,paste("mseweight_",this.vname,".",outform[o],sep=""))
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
      #------------------------------------------------------------------------------------#


      #----- Map. -------------------------------------------------------------------------#
      image.map( x              = region$lon
               , y              = region$lat
               , z              = weighted.stat
               , colour.palette = this.cscheme
               , levels         = this.levels
               , nlevels        = this.nlevels
               , na.col         = miss.colour
               , xlim           = longitude.limit
               , ylim           = latitude.limit
               , x.axis.options = list(side=1,at=lonplot$at,labels=lonplot$labels)
               , y.axis.options = list(side=2,at=latplot$at,labels=latplot$labels,las=1)
               , main.title     = list( text = paste(this.desc," - Composite")
                                      , font = 2
                                      , cex  = cex.main
                                      )#end list
               , key.log        = FALSE
               , key.title      = title( main     = this.key.main
                                       , cex.main = 0.8 * cex.main
                                       , line     = 1
                                       )#end title
               , plot.after     = list( abline     = list( v   = lonplot$at
                                                         , h   = latplot$at
                                                         , lwd = 1
                                                         , lty = "dotted"
                                                         , col = grid.colour
                                                         )#end list
                                      , southammap = list(lwd=1,col=grey.fg   )
                                      , amazonmap  = list(lwd=2,col=foreground)
                                      , text       = list( x      = site$lon
                                                         , y      = site$lat
                                                         , labels = site$iata
                                                         , col    = site[[col.key]]
                                                         , cex    = 1.5*site$cex
                                                         , font   = site$fnt
                                                         )#end list
                                      )#end list
               , f.key          = 1/6
               )#end image.map
      #------------------------------------------------------------------------------------#



      #----- Close the device. ------------------------------------------------------------#
      if (outform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
      bye = clean.tmp()
      #------------------------------------------------------------------------------------#
   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(nstatvar))
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
cat (" + Plot the soil texture...","\n")
#------------------------------------------------------------------------------------------#
#      Find the plot annotations.                                                          #
#------------------------------------------------------------------------------------------#
letitre       = "Soil classes - ED.2.2"
lesetiquettes = c("Sand","Clay","Silt")
lalegende     = paste(soil.ed21$ntext,soil.ed21$name,sep=" - ")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Create the tags for the axes.                                                       #
#------------------------------------------------------------------------------------------#
this      = seq(from=0.2,to=0.8,by=0.2)
neg       = rep(x=-0.04,times=length(this))
fill      = 1. - this - neg
sand.at    = cbind( x = this, y = neg , z = fill )
clay.at    = cbind( x = fill, y = this, z = neg  )
silt.at    = cbind( x = neg , y = fill, z = this )
textlab    = seq(from=20,to=80,by=20)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Plot the soil texture map.                                                           #
#------------------------------------------------------------------------------------------#
for (o in sequence(nout)){
   #------ Open file. ---------------------------------------------------------------------#
   fichier = file.path(outroot,paste("amzmap_stext.",outform[o],sep=""))
   if(outform[o] == "x11"){
      X11(width=ssize$width,height=ssize$height,pointsize=ptsz)
   }else if(outform[o] == "png"){
      png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
         ,pointsize=ptsz,res=depth)
   }else if(outform[o] == "eps"){
      postscript(file=fichier,width=ssize$width,height=ssize$height
                ,pointsize=ptsz,paper=ssize$paper)
   }else if(outform[o] == "pdf"){
      pdf(file=fichier,onefile=FALSE
         ,width=ssize$width,height=ssize$height,pointsize=ptsz,paper=ssize$paper)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Split the data plotting region into 2.  We put the legend outside. --------------#
   par(par.user)
   layout(cbind(2,1),widths=c(7,3))
   #---------------------------------------------------------------------------------------#




   #----- Plot legend. --------------------------------------------------------------------#
   par(mar=c(2.1,2.1,2.1,2.1))
   triplot( x            = soil.ed21$xsand
          , y            = soil.ed21$xclay
          , z            = soil.ed21$xsilt
          , label        = c("","","")
          , frame        = TRUE
          , grid         = FALSE
          , pch          = 16
          , col          = soil.ed21$colour
          , cex          = 0.1
          )#end triplot
   #---------------------------------------------------------------------------------------#


   #----- Plot the axis names. ------------------------------------------------------------#
   triframe(label = lesetiquettes,label.col=foreground)
   #---------------------------------------------------------------------------------------#

   #----- Plot the axis labels. -----------------------------------------------------------#
   text(x=tritrafo(sand.at),labels=textlab,srt=-60,cex=0.7,col=foreground,adj=c(0.5,0.8))
   text(x=tritrafo(clay.at),labels=textlab,srt=  0,cex=0.7,col=foreground,adj=c(0.8,0.5))
   text(x=tritrafo(silt.at),labels=textlab,srt= 60,cex=0.7,col=foreground,adj=c(0.3,0.5))
   #---------------------------------------------------------------------------------------#



   #----- Plot the background colours as polygons. ----------------------------------------#
   for (n in sequence(n.soil.ed21)){
      nst = soil.ed21$ntext[n]
      if (any(region$stext == nst,na.rm=TRUE) || any(stext.pret == nst,na.rm=TRUE)){
         polygon( x   = tritrafo( x = stext.polygon[[nst]]$sand
                                , y = stext.polygon[[nst]]$clay
                                , z = stext.polygon[[nst]]$silt
                                )#end tritrafo
                , density = -1
                , col     = stext.cols[nst]
                , border  = "transparent"
                )#end polygon
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#

   #----- Plot the grid. ------------------------------------------------------------------#
   #trigrid(x=seq(from=0.05,to=0.95,by=0.05),lty="dashed",col=grid.colour)
   #---------------------------------------------------------------------------------------#


   #----- Plot the class edges. -----------------------------------------------------------#
   for (n in 1:nstext.lines){
      trilines( x   = stext.lines[[n]]$sand
              , y   = stext.lines[[n]]$clay
              , z   = stext.lines[[n]]$silt
              , col = grey.mg
              ,lwd=2)
   }#end for
   #---------------------------------------------------------------------------------------#


   #----- Plot the sites. -----------------------------------------------------------------#
   for (n in sequence(n.soil.ed21)){
      nst = soil.ed21$ntext[n]
      if (any(region$stext == nst,na.rm=TRUE) || any(stext.pret == nst,na.rm=TRUE)){
         text( x      = tritrafo( x = soil.ed21$xsand[n]
                                , y = soil.ed21$xclay[n]
                                , z = soil.ed21$xsilt[n]
                                )#end tritrafo
             , labels = soil.ed21$key[n]
             , col    = foreground
             , srt    = 0
             , cex    = 0.7*cex.ptsz
             , font   = 1
             )#end text
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Plot the map.                                                                     #
   #---------------------------------------------------------------------------------------#
   par(mar=c(3.1,3.1,4.1,0.1))
   plot.new()
   plot.window(xlim=longitude.limit,ylim=latitude.limit)
   axis(side=1,at=lonplot$at,labels=lonplot$labels)
   axis(side=2,at=latplot$at,labels=latplot$labels,las=1)
   rect(xleft=region$wlon,xright=region$elon,ybottom=region$slat,ytop=region$nlat
       ,col=stext.cols[region$stext],border=stext.cols[region$stext])
   abline(v=lonplot$at,h=latplot$at,lwd=1,lty="dotted",col=grid.colour)
   southammap(lwd=1,col=grey.fg)
   amazonmap(lwd=2,col=foreground)
   text(x=site$lon,y=site$lat,labels=site$iata,col=site$col1,cex=site$cex,font=site$fnt)
   text(x=ptref$lon,y=ptref$lat,labels=ptref$lab,col=ptref$col,cex=ptref$cex*cex.ptsz
       ,font=ptref$fnt,adj=c(0.5,0.5))
   title (main=letitre)
   box()
   #---------------------------------------------------------------------------------------#

   #----- Close the device. ---------------------------------------------------------------#
   if (outform[o] == "x11"){
      locator(n=1)
      dev.off()
   }else{
      dev.off()
   }#end if
   bye = clean.tmp()
   #---------------------------------------------------------------------------------------#
}#end for
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Plot the 3-D plot of temperature, rainfall, and ET, coloured by some plant community #
# variable.                                                                                #
#------------------------------------------------------------------------------------------#
cat (" + Plot the 3-D plot of temperature, rainfall, and ET...","\n")

for (v in sequence(nthree.dim)){
   #----- Get the info for this colour plot. ----------------------------------------------#
   vname.fourth = three.dim[[v]]$vname
   key.fourth   = three.dim[[v]]$key
   unit.fourth  = untab[[three.dim[[v]]$unit]]
   leg.fourth   = desc.unit(desc=key.fourth,unit=unit.fourth)
   cat ("   - ",key.fourth,"...","\n")
   #---------------------------------------------------------------------------------------#


   #----- Copy variables to some generic place holders. -----------------------------------#
   first.var  = c(ymean$atm.temp)
   #first.var  = c(ymean$rshort  )
   second.var = c(ymean$rain    )
   third.var  = c(ymean$et      )
   #---------------------------------------------------------------------------------------#



   #------ Axes labels. -------------------------------------------------------------------#
   lex   = desc.unit(desc="Air temperature"    ,unit=untab$degC )
   #lex   = desc.unit(desc="Shortwave radiation",unit=untab$wom2 )
   ley   = desc.unit(desc="Rainfall"           ,unit=untab$mmoyr)
   lez   = desc.unit(desc="Evapotranspiration" ,unit=untab$mmoyr)
   #---------------------------------------------------------------------------------------#


   #----- Split the data into colours. ----------------------------------------------------#
   pretty.fourth   = pretty(ymean[[vname.fourth]],n=n.colourbar)
   n.pretty.fourth = length(pretty.fourth)
   fourth.breaks   = c(-Inf,pretty.fourth[-c(1,n.pretty.fourth)],Inf)
   fourth.colour   = cut(ymean[[vname.fourth]],fourth.breaks)
   fourth.cscheme  = clife(n=length(levels(fourth.colour))-1)
   fourth.colour   = fourth.cscheme[match(fourth.colour,levels(fourth.colour))]
   #---------------------------------------------------------------------------------------#




   #----- Find the axes annotations and limits. -------------------------------------------#
   x.at     = pretty(first.var ,n=6)
   y.at     = pretty(second.var,n=6)
   z.at     = pretty(third.var ,n=6)
   x.labels = x.at
   y.labels = y.at
   z.labels = z.at
   nx.at    = length(x.at)
   ny.at    = length(y.at)
   nz.at    = length(z.at)
   xlimit   = range(c(first.var ))
   ylimit   = range(c(second.var))
   zlimit   = range(c(third.var ))
   xfloor   = sort(unique(c(x.at,xlimit)))
   yfloor   = sort(unique(c(y.at,ylimit)))
   zfloor   = 0 * outer(xfloor,yfloor) + min(c(z.at,zlimit))
   zceiling = 0 * outer(xfloor,yfloor) + max(c(z.at,zlimit))
   xlimit   = range(xfloor)
   ylimit   = range(yfloor)
   zlimit   = range(c(zfloor,zceiling,z.at))
   x.mid    = mean(xlimit)
   y.mid    = mean(ylimit)
   z.mid    = mean(zlimit)
   #---------------------------------------------------------------------------------------#



   for (s in sequence(n.stext.pret)){
      #----- Concatenate similar soils to the same plot. ----------------------------------#
      sel             = use.stext[c(ymean$stext)] == stext.pret
      sel[is.na(sel)] = FALSE
      letitre         = paste("Annual means",stext.names[stext.pret[s]],sep=" - ")
      suffix          = paste("stext",sprintf("%2.2i",stext.pret[s]),sep="")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Don't enter here unless there is at least one point to plot...                #
      #------------------------------------------------------------------------------------#
      if (any(sel)){
         #----- Get only the data for this soil type. -------------------------------------#
         sel        = sample(which(sel))
         first.now  = first.var     [sel]
         second.now = second.var    [sel]
         third.now  = third.var     [sel]
         fourth.now = fourth.colour [sel]
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Loop over output formats.                                                   #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #------ Open file. ------------------------------------------------------------#
            fichier = file.path(outroot,paste("ettrplot_",vname.fourth,"_",suffix
                               ,".",outform[o],sep=""))
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
            #------------------------------------------------------------------------------#



            #----- Split the data plotting region into 2.  We put the legend outside. -----#
            par(par.user)
            layout(cbind(2,1),widths=c(6,1))
            #------------------------------------------------------------------------------#




            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(4.1,0.5,4.1,3.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=range(pretty.fourth),xaxs="i",yaxs="i")
            rect( xleft   = 0
                , ybottom = pretty.fourth[-n.pretty.fourth]
                , xright  = 1
                , ytop    = pretty.fourth[-1]
                , col     = fourth.cscheme
                , border  = "transparent"
                )#end rect
            axis (side=4,las=1)
            mtext(text=leg.fourth,side=3,cex.main=.7,line=.5)
            box()
            #------------------------------------------------------------------------------#




            #----- Don't plot floating stuff, just the projections onto the walls. --------#
            par(mar=c(1.1,1.1,4.1,1.1))
            pout = perspx( x        = xfloor
                         , y        = yfloor
                         , z        = zfloor
                         , xlim     = xlimit
                         , ylim     = ylimit
                         , zlim     = zlimit
                         , theta    = theta
                         , phi      = phi
                         , col      = if(ibackground==0){"grey99"}else{"grey01"}
                         , expand   = expz
                         , ticktype = "detailed"
                         , border   = NA
                         , shade    = shade
                         , ltheta   = ltheta
                         , axes     = FALSE
                         )#end perspx
            xxx  =  acos(cos(phi*pio180)*cos(theta*pio180)) / pio180
            yyy  =  xxx - 90
            zzz  =  - 90 + asin(expz*sin(phi*pio180)*cos(theta*pio180)) / pio180
            mtext(side=3,text=letitre,line=0.5,font=2,cex=cex.main)
            #------------------------------------------------------------------------------#



            #----- Add axes. --------------------------------------------------------------#
            paxis3d(edge="X--",pmat=pout,at=x.at,cex=.9*cex.ptsz,labels=x.labels,labdist=.2)
            paxis3d(edge="Y--",pmat=pout,at=y.at,cex=.9*cex.ptsz,labels=y.labels,labdist=.2)
            paxis3d(edge="Z-+",pmat=pout,at=z.at,cex=.9*cex.ptsz,labels=z.labels,labdist=.2)
            mtext3d(edge="X--",pmat=pout,labels=lex,cex=cex.ptsz,srt=xxx,at=x.mid,dist=.4)
            mtext3d(edge="Y--",pmat=pout,labels=ley,cex=cex.ptsz,srt=yyy,at=y.mid,dist=.4)
            mtext3d(edge="Z-+",pmat=pout,labels=lez,cex=cex.ptsz,srt=zzz,at=z.mid,dist=.4)
            #---------------------------------------------------------------------------------#



            #----- Grid lines. ------------------------------------------------------------#
            x.at.3 = rep(x.at,each=3)
            y.at.3 = rep(y.at,each=3)
            z.at.3 = rep(z.at,each=3)
            xl.3   = rep(xlimit[2],each=3)
            yl.3   = rep(ylimit[2],each=3)
            zl.3   = rep(zlimit[1],each=3)

            glist = list( trans3d(x=rep(c(xlimit,NA),times=ny.at),y=y.at.3,z=zl.3,pmat=pout)
                        , trans3d(x=rep(c(xlimit,NA),times=nz.at),y=yl.3,z=z.at.3,pmat=pout)
                        , trans3d(x=x.at.3,y=rep(c(ylimit,NA),times=nx.at),z=zl.3,pmat=pout)
                        , trans3d(x=xl.3,y=rep(c(ylimit,NA),times=nz.at),z=z.at.3,pmat=pout)
                        , trans3d(x=x.at.3,y=yl.3,z=rep(c(zlimit,NA),times=nx.at),pmat=pout)
                        , trans3d(x=xl.3,y=y.at.3,z=rep(c(zlimit,NA),times=ny.at),pmat=pout)
                        )#end list
            for (gg in sequence(length(glist))){
               lines(x = glist[[gg]],lty="dotted",col=grid.colour,lwd=1.5)
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Decide whether to plot the projections to the walls or the the floating  #
            # points.                                                                      #
            #------------------------------------------------------------------------------#
            if (proj.wall){

               #---------------------------------------------------------------------------#
               #      Plot line projection on the X = max(X) plane.                        #
               #---------------------------------------------------------------------------#
               points( trans3d( x    = 0. * first.now + xlimit[2]
                              , y    = second.now
                              , z    = third.now
                              , pmat = pout 
                              )#end trans3d
                     , type = "p"
                     , pch  = 16
                     , cex  = cex.pt3d
                     , col  = fourth.now
                     )#end lines
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Plot line projection on the Y = max(Y) plane.                        #
               #---------------------------------------------------------------------------#
               points( trans3d( x    = first.now
                              , y    = 0. * second.now + ylimit[2]
                              , z    = third.now
                              , pmat = pout 
                              )#end trans3d
                     , type = "p"
                     , pch  = 16
                     , cex  = cex.pt3d
                     , col  = fourth.now
                     )#end lines
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Plot line projection on the Z = min(Z) plane.                        #
               #---------------------------------------------------------------------------#
               points( trans3d( x    = first.now
                              , y    = second.now
                              , z    = 0 * third.now + zlimit[1]
                              , pmat = pout 
                              )#end trans3d
                     , type = "p"
                     , pch  = 16
                     , cex  = cex.pt3d
                     , col  = fourth.now
                     )#end lines
               #---------------------------------------------------------------------------#
            }else{
               #---------------------------------------------------------------------------#
               #      Plot floating objects.                                               #
               #---------------------------------------------------------------------------#
               points( trans3d( x    = first.now
                              , y    = second.now
                              , z    = third.now
                              , pmat = pout 
                              )#end trans3d
                     , type = "p"
                     , pch  = 16
                     , cex  = cex.pt3d
                     , col  = fourth.now
                     )#end lines
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            bye = clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (any(sel))
      #------------------------------------------------------------------------------------#
   }#end for (s in sequence(n.stext.pret))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(nthree.dim))
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Plot the box plot of ET0/Rain ratio by drought length for each soil texture.         #
#------------------------------------------------------------------------------------------#
cat (" + Plot the box plot of the ET0/Rain ratio...","\n")
for (s in sequence(n.stext.pret)){
   #----- Concatenate similar soils to the same plot. -------------------------------------#
   sel             = use.stext[c(ymean$stext)] == stext.pret
   sel[is.na(sel)] = FALSE
   letitre         = stext.names[stext.pret[s]]
   suffix          = paste("stext",sprintf("%2.2i",stext.pret[s]),sep="")
   #---------------------------------------------------------------------------------------#


   #----- Select this soil texture. -------------------------------------------------------#
   this        = drought[[suffix]]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether there is at least one point to be plotted.                          #
   #---------------------------------------------------------------------------------------#
   if ( length(this$dlen) > 0){

      dlen.breaks = c(seq(from=0.5,to=12.5,by=1),seq(from=17.5,to=37.5,by=5),Inf)
      dlen.lab    = 0.5 * ( dlen.breaks[-1] + dlen.breaks[-length(dlen.breaks)])
      dlen.lab[length(dlen.lab)] = dlen.lab[length(dlen.lab)-1] + 5.0
      idx         = cut(this$dlen,dlen.breaks)
      bp          = split(x=this$et0/this$rain,f=idx,drop=FALSE)
      names(bp)   = dlen.lab
      xlimit      = c(0,length(bp))+0.5
      xat         = sequence(length(bp))
      xgrid       = sequence(length(bp)+1) - 0.5
      ylimit      = range(unlist(bp),finite=TRUE)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Loop over output formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #------ Open file. ---------------------------------------------------------------#
         fichier = file.path(outroot,paste("boxplot_et0orain_",suffix
                            ,".",outform[o],sep=""))
         if(outform[o] == "x11"){
            X11(width=fsize$width,height=fsize$height,pointsize=ptsz)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=fsize$width,height=fsize$height
                      ,pointsize=ptsz,paper=fsize$paper)
         }else if(outform[o] == "pdf"){
            pdf(file=fichier,onefile=FALSE
               ,width=fsize$width,height=fsize$height,pointsize=ptsz,paper=fsize$paper)
         }#end if
         #---------------------------------------------------------------------------------#




         #----- Plot legend. --------------------------------------------------------------#
         par(par.user)
         plot.new()
         plot.window(xlim=xlimit,ylim=ylimit)
         abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="dotted")
         abline(h=1,col="chartreuse3",lwd=2,lty="dotdash")
         axis(side=1,at=xat,labels=dlen.lab)
         axis(side=2)
         title(main=letitre
              ,ylab=expression((dot(epsilon)[0]+dot(tau)[0])/r)
              ,xlab=desc.unit(desc="Drought length",unit=untab$month)
              )#end title
         boxplot(bp,col="slateblue3",add=TRUE,axes=FALSE)
         box()
         #---------------------------------------------------------------------------------#


         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         bye = clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
   }#end if (length(this$dlen) > 0)
   #---------------------------------------------------------------------------------------#
}#end for (s in sequence(n.stext.pret))
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Plot the skill plot for all precipitation data sets.                                #
#------------------------------------------------------------------------------------------#
cat (" + Plot the skill plot...","\n")
#----- Get the data for skill plot. -------------------------------------------------------#
use                = match(datrain$vname,dimnames(wmo$mod)[[3]])
obs.sdev           = array( data = apply(X=wmo$sdev,MARGIN=1,FUN=commonest,na.rm=TRUE)
                          , dim  = c(length(wmo$lon),ndatrain)
                          , dimnames = list(names(wmo$lon),datrain$vname)
                          )#end array
skill.obs          = wmo$obs
skill.mod          = wmo$mod[,,use]
skill.res          = wmo$res[,,use]
#----- Count comparable years, and discard stations that have very few measurements. ------#
skill.n            = apply(X=is.finite(wmo$res[,,use]),MARGIN=c(1,3),FUN=sum)
skill.n            = aperm( a    = array( data = skill.n
                                        , dim  = c(dim(skill.n),dim(skill.mod)[2])
                                        )#end array
                          , perm = c(1,3,2)
                          )#end aperm
dimnames(skill.n)  = dimnames(skill.mod)
discard            = skill.n < n.lsq.min
skill.mod[discard] = NA
skill.res[discard] = NA
skill.n  [discard] = NA
#------------------------------------------------------------------------------------------#




#----- Standardise the bias and scattering of residuals. ----------------------------------#
mod.bias  = apply(X=skill.res,MARGIN=c(1,3),FUN=mean,na.rm=TRUE) / obs.sdev
mod.sigma = apply(X=skill.res,MARGIN=c(1,3),FUN=sd  ,na.rm=TRUE) / obs.sdev
mod.n     = apply(X=skill.n  ,MARGIN=c(1,3),FUN=commonest ,na.rm=TRUE)
mod.bias [!is.finite(mod.bias )] = NA
mod.sigma[!is.finite(mod.sigma)] = NA
#------------------------------------------------------------------------------------------#


#----- Retrieve bias and sigma to make the RMSE hemisphere. -------------------------------#
xy.range     = 1.04 * max(abs(c(mod.bias,mod.sigma)),na.rm=TRUE)
bias.range   = xy.range * c(-1,1)
sigma.range  = xy.range * c( 0,1)
r2.range     = range(1-xy.range^2,1)
#------------------------------------------------------------------------------------------#



#----- Retrieve bias and sigma to make the RMSE hemisphere. -------------------------------#
nyears.min  = min(c(skill.n),na.rm=TRUE)
nyears.max  = max(c(skill.n),na.rm=TRUE)
nyears.leg  = sort(unique(c(nyears.min,nyears.max,pretty(skill.n,n=4))))
keep        = nyears.leg %in% seq(from=nyears.min,to=nyears.max,by=1)
nyears.leg  = nyears.leg[keep]
skill.cex   = 0.5 + 1.5 * ( skill.n    - nyears.min ) / ( nyears.max - nyears.min )
skill.cex   = apply(X=skill.cex,MARGIN=c(1,3),FUN=commonest,na.rm=TRUE)
cex.leg     = 0.5 + 1.5 * ( nyears.leg - nyears.min ) / ( nyears.max - nyears.min )
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Loop over output formats.                                                            #
#------------------------------------------------------------------------------------------#
for (o in sequence(nout)){
   #------ Open file. ---------------------------------------------------------------------#
   fichier = file.path(outroot,paste("skill_rain.",outform[o],sep=""))
   if(outform[o] == "x11"){
      X11(width=fsize$width,height=fsize$height,pointsize=ptsz)
   }else if(outform[o] == "png"){
      png(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
         ,pointsize=ptsz,res=depth)
   }else if(outform[o] == "eps"){
      postscript(file=fichier,width=fsize$width,height=fsize$height
                ,pointsize=ptsz,paper=fsize$paper)
   }else if(outform[o] == "pdf"){
      pdf(file=fichier,onefile=FALSE
         ,width=fsize$width,height=fsize$height,pointsize=ptsz,paper=fsize$paper)
   }#end if
   #---------------------------------------------------------------------------------------#




   #----- Split plotting area into 2. -----------------------------------------------------#
   par(par.user)
   layout(mat=rbind(c(3,3),c(1,2)),heights=c(5,1),widths=c(2,1))
   #---------------------------------------------------------------------------------------#




   #----- Plot legend. --------------------------------------------------------------------#
   par(mar=c(0.1,0.1,0.1,0.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend( x      = "center"
         , inset  = 0.0
         , legend = datrain$desc
         , col    = datrain$colour
         , pch    = pch.skill
         , pt.lwd = 2
         , title  = expression(bold("Rainfall data"))
         , ncol   = min(3,pretty.box(ndatrain)$ncol)
         , cex    = 0.9 * cex.ptsz
         , xpd    = TRUE
         )#end legend
   #---------------------------------------------------------------------------------------#




   #----- Plot legend. --------------------------------------------------------------------#
   par(mar=c(0.1,0.1,0.1,0.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend( x      = "center"
         , inset  = 0.0
         , legend = nyears.leg
         , col    = foreground
         , pch    = pch.skill
         , pt.cex = cex.leg
         , pt.lwd = 2
         , title  = expression(bold("Years"))
         , ncol   = min(3,pretty.box(length(nyears.leg))$ncol)
         , cex    = 0.9 * cex.ptsz
         , xpd    = TRUE
         )#end legend
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Plot the skew plot.                                                              #
   #---------------------------------------------------------------------------------------#
   par(mar=c(5.1,4.1,4.1,2.1))
   myskill = NULL
   for (d in datrain$order){
      for (w in seq_along(wmo$lon)){
         #----- Current data set information. ---------------------------------------------#
         now.obs = skill.obs[w,  ]
         now.mod = skill.mod[w,,d]
         now.cex = skill.cex[w, d]
         now.col = datrain$colour[d]
         #---------------------------------------------------------------------------------#



         #----- Plot the skill. -----------------------------------------------------------#
         myskill = skill.plot( obs            = now.obs
                             , obs.options    = list(col=foreground,cex=2.0)
                             , mod            = now.mod
                             , mod.options    = list( col = now.col
                                                    , pch = pch.skill
                                                    , cex = now.cex
                                                    , lwd = 3
                                                    )#end list
                             , main           = "Annual rainfall"
                             , bias.lim       = bias.range
                             , r2.lim         = r2.range
                             , r2.options     = list(col=grid.colour)
                             , nobias.options = list(col=khaki.mg)
                             , rmse.options   = list( col = orange.mg
                                                    , lty = "dotdash"
                                                    , lwd = 2
                                                    , bg  = background)
                             , cex.xyzlab     = 1.4
                             , cex.xyzat      = 1.4
                             , skill          = myskill
                             , normalise      = TRUE
                             , mar            = c(5,4,4,3)+0.1
                             )#end skill.plot
         #---------------------------------------------------------------------------------#
      }#end for (d in seq_along(wmo$lon))
      #------------------------------------------------------------------------------------#
   }#end for (d in sequence(ndatrain))
   #---------------------------------------------------------------------------------------#


   #----- Close the device. ---------------------------------------------------------------#
   if (outform[o] == "x11"){
      locator(n=1)
      dev.off()
   }else{
      dev.off()
   }#end if
   bye = clean.tmp()
   #---------------------------------------------------------------------------------------#
}#end for (o in sequence(nout))
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Plot the number of years of each site used.                                         #
#------------------------------------------------------------------------------------------#
cat (" + Plot the WMO statistics scaled by data availability...","\n")

#----- Get the number of available years and find the size scale. -------------------------#
nyears.min  = min(wmo$nyears,na.rm=TRUE)
nyears.max  = max(wmo$nyears,na.rm=TRUE)
nyears.leg  = sort(unique(c(nyears.min,nyears.max,pretty(wmo$nyears,n=4))))
keep        = nyears.leg %in% seq(from=nyears.min,to=nyears.max,by=1)
nyears.leg  = nyears.leg[keep]
wmo$cex     = 0.7 + 1.3 * ( wmo$nyears - nyears.min ) / ( nyears.max - nyears.min )
cex.leg     = 0.7 + 1.3 * ( nyears.leg - nyears.min ) / ( nyears.max - nyears.min )
#------------------------------------------------------------------------------------------#


#----- Pack the legend options to a list. -------------------------------------------------#
xyz.legend  = list( x      = "center"
                  , inset  = 0.0
                  , legend = nyears.leg
                  , col    = foreground
                  , pch    = pch.wmo
                  , pt.cex = cex.leg
                  , title  = expression(bold("Available years"))
                  , ncol   = length(nyears.leg)
                  , cex    = 0.9 * cex.ptsz
                  , xpd    = TRUE
                  )#end legend
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
for (v in sequence(nstatvar)){
   this.stat     = statvar[[v]]
   this.vname    = this.stat$vname
   this.desc     = this.stat$desc
   this.key.main = desc.unit(desc=NULL,unit=untab[[this.stat$unit]])
   this.cscheme  = get(this.stat$cscheme)
   this.min      = this.stat$min
   this.max      = this.stat$max
   this.plog     = this.stat$plog
   plotit        = this.vname %in% names(wmo)

   if (plotit){

      #---- Select colour for sites. ------------------------------------------------------#
      if (this.stat$cscheme %in% col1.schemes){
         col.key = "col1"
      }else{
         col.key = "col2"
      }#end if
      #------------------------------------------------------------------------------------#



      cat ("   - ",this.desc,"...","\n")


      #---- Grab variable. ----------------------------------------------------------------#
      var.now = wmo[[this.vname]]
      #------------------------------------------------------------------------------------#



      #------ Limits for the z key scale. -------------------------------------------------#
      if (! ( is.null(this.min) || is.null(this.max) )){
         var.now      = pmax(this.min,pmin(this.max,var.now))
         this.limit   = c(this.min,this.max)
      }else{
         this.limit   = range(var.now,na.rm=TRUE)
      }#end if
      if (this.plog){
         this.levels  = pretty.log(this.limit,n=n.colourbar)
         this.levels  = c( min(this.levels) / mean(ediff(this.levels))
                         , this.levels
                         , max(this.levels) * mean(ediff(this.levels))
                         )#end c
         this.nlevels = length(this.levels)
      }else{
         this.levels  = pretty(this.limit,n=n.colourbar)
         this.levels  = c( min(this.levels) - mean(diff(this.levels))
                         , this.levels
                         , max(this.levels) + mean(diff(this.levels))
                         )#end c
         this.nlevels = length(this.levels)
      }#end if
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Loop over output formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #------ Open file. ---------------------------------------------------------------#
         fichier = file.path(outroot,paste("wmo_",this.vname,".",outform[o],sep=""))
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
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Plot the skew plot.                                                        #
         #---------------------------------------------------------------------------------#
         par(par.user)
         xyz.plot  ( x              = wmo$lon
                   , y              = wmo$lat
                   , z              = var.now
                   , xlim           = longitude.limit
                   , ylim           = latitude.limit
                   , colour.palette = this.cscheme
                   , levels         = this.levels
                   , nlevels        = this.nlevels
                   , pch            = pch.wmo
                   , cex            = wmo$cex
                   , xyz.before     = list( southammap = list(col=grey.mg)
                                          , amazonmap  = list(col=foreground,lwd=2.0)
                                          )#end list
                   , xyz.after      = list( text       = list( x      = site$lon
                                                             , y      = site$lat
                                                             , labels = site$iata
                                                             , col    = site[[col.key]]
                                                             , cex    = 1.5*site$cex
                                                             , font   = site$fnt
                                                             )#end list
                                          )#end list
                   , xyz.title      = list ( main = this.desc    )
                   , key.title      = title( main = this.key.main)
                   , xyz.xaxis      = list(side=1,at=lonplot$at,labels=lonplot$labels)
                   , xyz.yaxis      = list(side=2,at=latplot$at,labels=latplot$labels,las=1)
                   , xyz.legend     = xyz.legend
                   , key.width      = 5.0
                   , leg.height     = 8.0
                   )#end skill.plot
         #---------------------------------------------------------------------------------#


         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         bye = clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
   }#end if (plotit)
   #---------------------------------------------------------------------------------------#
}#end for(v in sequence(nstatvar))
#==========================================================================================#
#==========================================================================================#
