#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
here        = getwd()                                       #   Current directory
srcdir      = c( "/home/b14275/Util/Rsc"                    # Possible paths with libraries
               , "/Users/mlongo/Util/Rsc"                   #    R will select the first
               , "/prj/prjidfca/marcosl/Util/Rsc"           #    one that is found.
               , "/prj/bramsolam/marcos.longo/Util/Rsc"     #
               , "/scratch/bramsolam/marcos.longo/Util/Rsc" #
               , "/n/home00/mlongo/util/Rsc"                #
               )#end c                                      #
ibackground = 0                                             # Sought background
                                                            # 0 -- white
                                                            # 1 -- black
                                                            # 2 -- dark grey
#----- Output directory -------------------------------------------------------------------#
outroot = file.path(here,paste0("patch_pdf_ibg_",sprintf("%2.2i",ibackground)))
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.hour  = c(FALSE,TRUE)[2]
rdata.path   = file.path(here,"RData_patches")
rdata.suffix = "patches_ed22.RData"
#------------------------------------------------------------------------------------------#



#----- Default settings. ------------------------------------------------------------------#
emean.yeara = 2008  # First year
emean.yearz = 2010  # Last year
eshow.yeara = 2009  # First year to show
eshow.yearz = 2009  # First year to show
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Quantiles that determine which patches to show.                                       #
#------------------------------------------------------------------------------------------#
qshow = c(0,0.25,0.50,0.75,1.00)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Vertical profile variables.                                                           #
#------------------------------------------------------------------------------------------#
hgt.cdf.brks = seq(from=0,to=45,by=1)
hgt.lyr.brks = seq(from=0,to=45,by=5)
#------------------------------------------------------------------------------------------#



#----- Map settings. ----------------------------------------------------------------------#
pat.min.area = 0.0025
n.min.area   = 20
dxy.gap      = 20
#------------------------------------------------------------------------------------------#


#----- Original number of patches: this is used to determine the weighted variance. -------#
npat.orig = 4203
#------------------------------------------------------------------------------------------#


#----- Confidence band. ----------------------------------------------------------------------#
cband         = c(-1,2*pnorm(q=1)-1)[2]
cband.leg     = c(test="Heterogeneous",ctrl="Homogeneous")
cband.col     = c(test="#9EEEEE"      ,ctrl="#7F7F7F"    )
cmean.col     = c(test="#36A6A6"      ,ctrl="#404040"    )
cmean.lty     = c(test="twodash"      ,ctrl="solid"      )
cband.angle   = c(test=-45            ,ctrl=+45          )
cband.density = c(test= 60            ,ctrl= 60          )
cband.lwd     = c(test=1/2            ,ctrl= 1/2         )
#------------------------------------------------------------------------------------------#


#----- Limits for cumulative LAI layers. --------------------------------------------------#
clai.lim      = c(0,0.5,3.0,Inf)
clai.top.lim  = list(zupr=clai.lim[1],zmid=clai.lim[2],zlwr=clai.lim[3])
clai.bot.lim  = list(zupr=clai.lim[2],zmid=clai.lim[3],zlwr=clai.lim[4])
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Site settings:                                                                       #
#                                                                                          #
# eort      -- first letter ("e" or "t")                                                   #
#                                                                                          #
# sites     -- a list of the sites that could be processed (see use.sites). Each list      #
#              element is also a list with the following information:                      #
#              iata      -- 3-letter code to refer to this site                            #
#              desc      -- Full name of the site (for titles).                            #
#              pch       -- site symbols                                                   #
#              col       -- background colour                                              #
#              fg        -- foreground colour                                              #
#              drya      -- beginning of the dry season (MM/DD) -- NA skips dry season     #
#              dryz      -- end of the dry season       (MM/DD) -- NA skips dry season     #
#                           dryz can be less than drya, in which case two rectangles are   #
#                           plotted.                                                       #
#                                                                                          #
# use.sites -- controls which sites are actually run.  It does different things depending  #
#              on the type.                                                                #
#              Logical   -- If TRUE, it uses all sites; FALSE doesn't make sense.          #
#              Character -- vector with the 3-letter code ("iata") of selected sites       #
#              Numeric   -- Site index.                                                    #
#                                                                                          #
#------------------------------------------------------------------------------------------#
eort       = "t"
n          = 0
sites      = list()
n          = n + 1
sites[[n]] = list( iata = "tnf"
                 , desc = "Tapajos National Forest"
                 , pch  =  5
                 , col  = "#A3CC52"
                 , fg   = "#4B6614"
                 , drya = "07/13"
                 , dryz = "11/21"
                 )#end list
use.sites  = TRUE
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Eddy flux comparisons.                                                               #
#------------------------------------------------------------------------------------------#
n = 0
compvar       = list()
n             = n + 1
compvar[[ n]] = list( vnam     = "zupr.gpp"
                    , desc     = "GPP - Upper canopy"
                    , unit     = "kgwom2oday"
                    , cscheme  = "atlas"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "zmid.gpp"
                    , desc     = "GPP - Mid-canopy"
                    , unit     = "kgwom2oday"
                    , cscheme  = "atlas"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "zlwr.gpp"
                    , desc     = "GPP - Lower canopy"
                    , unit     = "kgwom2oday"
                    , cscheme  = "atlas"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "gpp"
                    , desc     = "Gross primary productivity"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "atlas"
                    , qmean    = TRUE
                    , clprof   = "gpp"
                    , clcum    = "nplant"
                    , hgtprof  = TRUE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "npp"
                    , desc     = "Net primary productivity"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "atlas"
                    , qmean    = TRUE
                    , clprof   = "npp"
                    , clcum    = "nplant"
                    , hgtprof  = TRUE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "plant.resp"
                    , desc     = "Plant respiration"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "iatlas"
                    , qmean    = TRUE
                    , clprof   = "plant.resp"
                    , clcum    = "nplant"
                    , hgtprof  = TRUE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "cba"
                    , desc     = "Carbon balance"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "atlas"
                    , qmean    = FALSE
                    , clprof   = "cba"
                    , clcum    = NA_character_
                    , hgtprof  = TRUE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "reco"
                    , desc     = "Ecosystem respiration"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "iatlas"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "nep"
                    , desc     = "Net Ecosystem Productivity"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "atlas"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "hflxca"
                    , desc     = "Sensible heat flux"
                    , unit     = "wom2"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "wflxca"
                    , desc     = "Water vapour flux"
                    , unit     = "kgwom2oday"
                    , cscheme  = "ipanoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "transp"
                    , desc     = "Transpiration"
                    , unit     = "kgwom2oday"
                    , cscheme  = "ipanoply"
                    , qmean    = TRUE
                    , clprof   = "transp"
                    , clcum    = "nplant"
                    , hgtprof  = TRUE
                    , yrcumul  = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "rshortup"
                    , desc     = "Upward SW radiation"
                    , unit     = "wom2"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "rlongup"
                    , desc     = "Upward LW radiation"
                    , unit     = "wom2"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "parup"
                    , desc     = "Upward PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.gnd"
                    , desc     = "Ground absorption - PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "rshort.gnd"
                    , desc     = "Ground absorption - SW"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "sm.stress"
                    , desc     = "Soil moisture stress"
                    , unit     = "empty"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.gpp"
                    , desc     = "Leaf GPP"
                    , unit     = "kgcom2loyr"
                    , cscheme  = "atlas"
                    , qmean    = TRUE
                    , clprof   = "leaf.gpp"
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.temp"
                    , desc     = "Mean Leaf Temperature"
                    , unit     = "degC"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.temp"
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.vpd"
                    , desc     = "Mean Leaf VPD"
                    , unit     = "hpa"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.vpd"
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.gsw"
                    , desc     = "Stomatal conductance"
                    , unit     = "kgwom2loday"
                    , cscheme  = "ipanoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.gsw"
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.leaf"
                    , desc     = "Leaf Absorption - PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.par"
                    , clcum    = "lai"
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.leaf.beam"
                    , desc     = "Leaf Absorption - Direct PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.par.beam"
                    , clcum    = "lai"
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.leaf.diff"
                    , desc     = "Leaf Absorption - Diffuse PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.par.diff"
                    , clcum    = "lai"
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.par"
                    , desc     = "Norm. Leaf Absorption - PAR"
                    , unit     = "umolom2los"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.par"
                    , clcum    = NA_character_
                    , hgtprof  = TRUE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.par.beam"
                    , desc     = "Norm. Leaf Absorption - Direct PAR"
                    , unit     = "umolom2los"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.par.beam"
                    , clcum    = NA_character_
                    , hgtprof  = TRUE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.par.diff"
                    , desc     = "Norm. Leaf Absorption - Diffuse PAR"
                    , unit     = "umolom2los"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    , clprof   = "leaf.par.diff"
                    , clcum    = NA_character_
                    , hgtprof  = TRUE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "assim.light"
                    , desc     = "Light-limited Assimilation"
                    , unit     = "umolom2los"
                    , cscheme  = "atlas"
                    , qmean    = TRUE
                    , clprof   = "assim.light"
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "assim.rubp"
                    , desc     = "RuBP-limited Assimilation"
                    , unit     = "umolom2los"
                    , cscheme  = "atlas"
                    , qmean    = TRUE
                    , clprof   = "assim.rubp"
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.lpar"
                    , desc     = "Daytime PAR absorportion by leaves"
                    , unit     = "hpa"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.lvpd"
                    , desc     = "Daytime Leaf VPD"
                    , unit     = "hpa"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.ltemp"
                    , desc     = "Daytime Leaf Temperature"
                    , unit     = "degC"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.sms"
                    , desc     = "Daytime soil moisture stress"
                    , unit     = "empty"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.lgsw"
                    , desc     = "Daytime stomatal conductance"
                    , unit     = "kgwom2loday"
                    , cscheme  = "ipanoply"
                    , qmean    = FALSE
                    , clprof   = NA_character_
                    , clcum    = NA_character_
                    , hgtprof  = FALSE
                    , yrcumul  = FALSE
                    )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.suffix  = "imetrad04"
sim.struct  = c(ctrl = "ihrzrad04", test = "ihrzrad02")
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform           = c("pdf")             # Formats for output file.  Supported format:
                                         #   - "X11" - for printing on screen
                                         #   - "eps" - for postscript printing
                                         #   - "png" - for PNG printing
                                         #   - "tif" - for TIFF printing
                                         #   - "pdf" - for PDF printing

byeold            = TRUE                 # Remove old files of the given format?

depth             = 600                  # PNG resolution, in pixels per inch
paper             = "square"             # Paper size, to define the plot shape
mpaper            = "double"             # Paper size for maps, to define the plot shape
ptsz              = 17                   # Font size.
f.leg             = 1/6                  # Factor to expand plot devices
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Switch controls to plot only the needed ones.                                       #
#------------------------------------------------------------------------------------------#
plot.tseries    = c(FALSE,TRUE)[2]
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Number of levels for colour palette.                                                #
#------------------------------------------------------------------------------------------#
ncolpal = 200
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE.    #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#



#----- Load some packages. ----------------------------------------------------------------#
if (any(file.exists(srcdir))){
   srcdir = (srcdir[file.exists(srcdir)])[1]
   source(file.path(srcdir,"load.everything.r"))
}else{
   stop("None of the paths provided in variable \"srcdir\" exists.")
}#end if (any(file.exists(srcdir)))
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#

#----- Create directory with RData. -------------------------------------------------------#
if (! file.exists(rdata.path)) dir.create(rdata.path)
#------------------------------------------------------------------------------------------#



#----- Count and configure quantiles. -----------------------------------------------------#
nquant     = length(qshow)
qkey.show  = paste0("q",sprintf("%3.3i",100.*qshow))
qdesc.show = paste0("Quantile: ",100*qshow,"%")
#------------------------------------------------------------------------------------------#


#----- Number of breaks for colour palette. -----------------------------------------------#
ncpbks     = ncolpal + 1
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Vertical profile variables.                                                           #
#------------------------------------------------------------------------------------------#
hgt.cdf.range             = range(hgt.cdf.brks)
hgt.lyr.range             = range(hgt.lyr.brks)
hgt.cdf.levs              = mid.points(hgt.cdf.brks)
hgt.lyr.levs              = mid.points(hgt.lyr.brks)
hgt.cdf.labels            = paste0("z",sprintf("%05.1f",hgt.cdf.levs))
hgt.lyr.labels            = paste0("z",sprintf("%05.1f",hgt.lyr.levs))
n.hgt.cdf                 = length(hgt.cdf.levs)
n.hgt.lyr                 = length(hgt.lyr.levs)
dhgt.cdf                  = mean(diff(hgt.cdf.brks))
dhgt.lyr                  = mean(diff(hgt.lyr.brks))
hgt.cdf.brks[n.hgt.cdf+1] = Inf
hgt.lyr.brks[n.hgt.lyr+1] = Inf
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Convert sites to data frame, and keep only those to be run.                         #
#------------------------------------------------------------------------------------------#
sites = list.2.data.frame(sites)
if (is.logical(use.sites)){
   use.sites = which(rep(use.sites,times=nrow(sites))[sequence(nrow(sites))])
}else if (is.character(use.sites)){
   use.sites = match(use.sites,sites$iata)
   use.sites = use.sites[! is.na(use.sites)]
}#end if
if (length(use.sites) == 0){
   stop ("You must specify at least one valid site")
}else{
   sites = sites[use.sites,,drop=FALSE]
}#end if
nsites = nrow(sites)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      List the variables.                                                                 #
#------------------------------------------------------------------------------------------#
compvar  = list.2.data.frame(compvar)
ncompvar = nrow(compvar)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Light settings.                                                                     #
#------------------------------------------------------------------------------------------#
light.key = c("Bright"     ,"Intermediate","Dark"  ,"Average","Homogeneous")
light.col = c("deepskyblue","dodgerblue3" ,"grey16","gold"   ,"firebrick3" )
light.lty = c("solid"      ,"solid"       ,"solid" ,"twodash","dotdash"    )
light.lwd = c(2.0          ,2.0           ,2.0     ,2.8      ,2.8          )
nlight    = length(light.key)
nltype    = nlight - 2
nwl       = nlight - 1 
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (1. - f.leg)
size0   = plotsize(proje=FALSE,paper=paper)
height0 = size0$height
width0  = size0$width
xsize   = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (2. - f.leg)
msize   = plotsize(proje=FALSE,paper=mpaper,extendfc="lon",extfactor=f.ext/2)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Define some utility functions to determine derived patch variables.                 #
#------------------------------------------------------------------------------------------#
sum.soil.c.fun = function(fast,slow,struct) fast + slow + struct
bowen.fun      = function(h,qw) ifelse(test = qw %!=% 0., yes = h/qw, no = NA)
tratio.fun     = function(tp,w) ifelse(test =  w %!=% 0., yes = tp/w, no = NA)
discard.fun    = function(x) x * NA
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Function that integrates variables by LAI layers within each cohort.                #
#------------------------------------------------------------------------------------------#
layer.gpp.list.fun = function(nplant,area,gpp,lai,ipa,top,bot){
   #----- Create data frame with nplant, gpp, and LAI.
   dat        = data.frame(nplant = nplant, area = area, gpp = gpp, lai = lai)
   dlist      = split(x=dat,f=ipa)
   ans        = sapply(X=dlist,FUN=layer.gpp.one.fun,top=top,bot=bot)
   names(ans) = NULL
   return(ans)
}#end function layer.gpp.list.fun
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Function that integrates variables by LAI layers within each cohort.                #
#------------------------------------------------------------------------------------------#
layer.gpp.one.fun = function(dat,top,bot){

   #------ We can't process empty data frames. --------------------------------------------#
   if (nrow(dat) > 0){
      #----- Discard empty layers. --------------------------------------------------------#
      dat = dat[dat$lai %>% 0,,drop=FALSE]
      #------------------------------------------------------------------------------------#
   }#end if (nrow(dat) > 0)
   #---------------------------------------------------------------------------------------#


   #----- Check again if the number of rows is greater than zero. -------------------------#
   if (nrow(dat) > 0){

      #----- Copy data to local variables. ------------------------------------------------#
      nplant = dat$nplant / dat$area
      gpp    = dat$gpp
      lai    = dat$lai
      #------------------------------------------------------------------------------------#


      #----- Cumulative LAI. --------------------------------------------------------------#
      clai     = cumsum(c(0,lai))
      nlai     = length(clai)
      clai.top = clai[-nlai]
      clai.bot = clai[   -1]
      #------------------------------------------------------------------------------------#



      #----- Integrate the productivity. --------------------------------------------------#
      fmult = ifelse( test = clai.bot > top & clai.top < bot
                    , yes  = (pmin(clai.bot,bot)-pmax(clai.top,top))/(clai.bot-clai.top)
                    , no   = 0.
                    )#end ifelse
      #------------------------------------------------------------------------------------#

      #----- Scale cohorts with fmult and return the sum. ---------------------------------#
      ans = sum(nplant * gpp * fmult)
      #------------------------------------------------------------------------------------#
   }else{
      #------ Nothing to return. ----------------------------------------------------------#
      ans = 0
      #------------------------------------------------------------------------------------#
   }#end if (nrow(dat) > 0)
   #---------------------------------------------------------------------------------------#


   #----- Return the sum. -----------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end slice lai.fun
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Function that integrates variables by LAI layers within each cohort.                #
#------------------------------------------------------------------------------------------#
layer.gpp.one.fun = function(dat,top,bot){

   #------ We can't process empty data frames. --------------------------------------------#
   if (nrow(dat) > 0){
      #----- Discard empty layers. --------------------------------------------------------#
      dat = dat[dat$lai %>% 0,,drop=FALSE]
      #------------------------------------------------------------------------------------#
   }#end if (nrow(dat) > 0)
   #---------------------------------------------------------------------------------------#


   #----- Check again if the number of rows is greater than zero. -------------------------#
   if (nrow(dat) > 0){

      #----- Copy data to local variables. ------------------------------------------------#
      nplant = dat$nplant / dat$area
      gpp    = dat$gpp
      lai    = dat$lai
      #------------------------------------------------------------------------------------#


      #----- Cumulative LAI. --------------------------------------------------------------#
      clai     = cumsum(c(0,lai))
      nlai     = length(clai)
      clai.top = clai[-nlai]
      clai.bot = clai[   -1]
      #------------------------------------------------------------------------------------#



      #----- Integrate the productivity. --------------------------------------------------#
      fmult = ifelse( test = clai.bot > top & clai.top < bot
                    , yes  = (pmin(clai.bot,bot)-pmax(clai.top,top))/(clai.bot-clai.top)
                    , no   = 0.
                    )#end ifelse
      #------------------------------------------------------------------------------------#

      #----- Scale cohorts with fmult and return the sum. ---------------------------------#
      ans = sum(nplant * gpp * fmult)
      #------------------------------------------------------------------------------------#
   }else{
      #------ Nothing to return. ----------------------------------------------------------#
      ans = 0
      #------------------------------------------------------------------------------------#
   }#end if (nrow(dat) > 0)
   #---------------------------------------------------------------------------------------#


   #----- Return the sum. -----------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end slice lai.fun
#------------------------------------------------------------------------------------------#


#------ Find out whether to plot annual means (it must have at least 2 years). ------------#
plot.ymean = eshow.yearz > eshow.yeara
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
outsimul   = file.path(outroot,sim.suffix)
outymean   = file.path(outsimul,"tseries_ymean")
outemean   = file.path(outsimul,"tseries_emean")
outcmean   = file.path(outsimul,"tseries_cmean")
outqmean   = file.path(outsimul,"tseries_qmean")
outcumlai  = file.path(outsimul,"cumlai_prof"  )
outhgtprof = file.path(outsimul,"height_prof"  )
outhgtcum  = file.path(outsimul,"height_cumsum")
outmmean   = file.path(outsimul,"maps_emean"   )
if (! file.exists(outroot   )             ) dir.create(outroot   )
if (! file.exists(outsimul  )             ) dir.create(outsimul  )
if (! file.exists(outymean  ) & plot.ymean) dir.create(outymean  )
if (! file.exists(outemean  )             ) dir.create(outemean  )
if (! file.exists(outcmean  )             ) dir.create(outcmean  )
if (! file.exists(outqmean  )             ) dir.create(outqmean  )
if (! file.exists(outcumlai )             ) dir.create(outcumlai )
if (! file.exists(outhgtprof)             ) dir.create(outhgtprof)
if (! file.exists(outhgtcum )             ) dir.create(outhgtcum )
if (! file.exists(outmmean  )             ) dir.create(outmmean  )
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Find limits for final display.                                                      #
#------------------------------------------------------------------------------------------#
when.show = chron(paste(c(1,12),c(1,1),c(eshow.yeara,eshow.yearz),sep="/"))
#------------------------------------------------------------------------------------------#



#----- Map settings. ----------------------------------------------------------------------#
n.xy     = ceiling(sqrt(n.min.area / pat.min.area))
n.gaps   = n.xy*n.xy
xy.gaps  = (sequence(n.xy)-0.5)*dxy.gap
xy.range = c(0,n.xy+1)*dxy.gap
#------------------------------------------------------------------------------------------#




#----- Confidence band. -------------------------------------------------------------------#
if (cband > 0){
   qlwr   = 0.5 - 0.5 * cband
   qupr   = 0.5 + 0.5 * cband
   cbdesc = paste0(sprintf("%2i",round(100*cband)),"% range")
}else{
   cbdesc = paste0("+/-",sprintf("%.1f",abs(cband))," S.D.")
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Retrieve all data that already exists.                                              #
#------------------------------------------------------------------------------------------#
res        = list()
loop.sites = integer(length=0)
for (s in sequence(nsites)){
   iata = sites$iata[s]
   #----- Find file name. -----------------------------------------------------------------#
   rdata.iata = file.path(rdata.path,paste(iata,sim.suffix,rdata.suffix,sep="_"))
   if (reload.hour && file.exists(rdata.iata)){
      #----- Reload data and copy to the general list. ------------------------------------#
      cat0(" + Load data from ",basename(rdata.iata),".")
      dummy       = load(file=rdata.iata)
      res.iata    = paste("res",iata,sep=".")
      res[[iata]] = get(res.iata)
      rm(res.iata)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Add site to read list. -------------------------------------------------------#
      loop.sites = c(loop.sites,s)
      #------------------------------------------------------------------------------------#
   }#end if (reload.hour && file.exists(rdata.iata))
   #---------------------------------------------------------------------------------------#
}#end for (s in sequence(nsites))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Loop over all sites, variables, and simulations, to prepare the data for the        #
# model comparison.  We want to use only times with actual measurements, so we will        #
# discard model results from times with no observation so all derived quantities have      #
# the same number of defined points (so if measurements are biased towards daytime, the    #
# model will also be equally biased).                                                      #
#------------------------------------------------------------------------------------------#
if (length(loop.sites) != 0) cat0(" + Process missing hourly data.")
for (s in loop.sites){
   #----- Get the basic information. ------------------------------------------------------#
   iata          = sites$iata[s]
   im            = match(iata,poilist$iata)
   this          = list()
   this$short    = poilist$short   [im]
   this$longname = sites$desc      [ s]
   this$iata     = poilist$iata    [im]
   this$lon      = poilist$lon     [im]
   this$lat      = poilist$lat     [im]
   this.sim      = list()
   cat0("   - Site :",this$longname,".")
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Get all the statistics and actual values for every simulation.                    #
   #---------------------------------------------------------------------------------------#
   cat0("    * Aggregate and find statistics for simulations for this site.")

   #----- Load hourly averages (Control). -------------------------------------------------#
   ctrl.name = paste0(eort,iata,"_",sim.struct["ctrl"],"_",sim.suffix)
   cat0("      # Control: ",ctrl.name,".")
   ctrl.path = file.path(here,ctrl.name)
   ctrl.file = file.path(ctrl.path,"rdata_month",paste0(ctrl.name,".RData"))
   load(ctrl.file)
   ctrl      = datum
   rm(datum)
   #---------------------------------------------------------------------------------------#


   #----- Load hourly averages. -----------------------------------------------------------#
   test.name = paste0(eort,iata,"_",sim.struct["test"],"_",sim.suffix)
   cat0("      # Test: ",test.name,".")
   test.path = file.path(here,test.name)
   test.file = file.path(test.path,"rdata_month",paste0(test.name,".RData"))
   load(test.file)
   test      = datum
   rm(datum)
   #---------------------------------------------------------------------------------------#


   #----- Read in the patch table. --------------------------------------------------------#
   base.ptable     = paste0("pix_ptable_isi001_",emean.yeara,"-01.txt")
   test.ptable     = file.path( test.path,"shade",base.ptable)
   test.ptable.bz2 = paste0(test.ptable,".bz2")
   test.ptable.gz  = paste0(test.ptable,".gz" )
   if (file.exists(test.ptable)){
      cat0("      # Test table: ",basename(test.name),".")
      ptable       = read.table(test.ptable,header=TRUE)
   }else if (file.exists(test.ptable.bz2)){
      temp.ptable  = file.path(tempdir(),base.ptable)
      dummy        = bunzip2(filename=test.ptable.bz2,destname=temp.ptable,remove=FALSE)
      ptable       = read.table(temp.ptable,header=TRUE)
      dummy        = file.remove(temp.ptable)
   }else if (file.exists(test.ptable.gz )){
      temp.ptable  = file.path(tempdir(),base.ptable)
      dummy        = gunzip(filename=test.ptable.gz,destname=temp.ptable,remove=FALSE)
      ptable       = read.table(temp.ptable,header=TRUE)
      dummy        = file.remove(temp.ptable)
   }else{
      stop(paste0(" File ",base.ptable," (or its gz/bz2 version) wasn't found!!!"))
   }#end if
   names(ptable) = tolower(gsub(pattern="_",replacement=".",x=names(ptable)))
   #---------------------------------------------------------------------------------------#




   #----- Create some variables to describe season and time of the day. -------------------#
   model = list()
   #----- Create time stamp for annual and monthly means. ---------------------------------#
   nyear         = emean.yearz - emean.yeara + 1
   model$toyear  = rep(seq(from=emean.yeara,to=emean.yearz,by=1),each=12)
   model$tomonth = chron( paste( month = rep( sequence(12), times=nyear   )
                               , day   = rep(            1, times=nyear*12)
                               , year  = model$toyear
                               , sep   = "/"
                               )#end paste
                        )#end chron
   nemean        = length(model$tomonth)
   model$imap    = sample.int( n       = nrow(ptable)
                             , size    = n.gaps
                             , replace = TRUE
                             , prob    = ptable$cci.area
                             )#end sample.int
   #---------------------------------------------------------------------------------------#





   #---- Save selection indices. ----------------------------------------------------------#
   sel        = ctrl$when %wr% model$tomonth
   eloop.ctrl = match(ctrl$when[sel],model$tomonth)
   sel        = test$when %wr% model$tomonth
   eloop.test = match(test$when[sel],model$tomonth)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Find number of mean diurnal cycle.                                               #
   #---------------------------------------------------------------------------------------#
   ndcycle       = ctrl$ndcycle
   dcstep        = ndcycle / day.hr
   dchour        = seq(from=0,by=dcstep,length.out=ndcycle)
   dclabel       = paste0(sprintf("%2.2i",dchour),"Z")
   model$ndcycle = ndcycle
   model$dchour  = dchour
   model$dclabel = dclabel
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Find some variables that must be always loaded.                                  #
   #---------------------------------------------------------------------------------------#
   npatches = length(ctrl$patch$age[[1]])
   if (length(unique(ptable$ipa)) != npatches){
      stop("Control and test must have the same number of patches.")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Define some patch characteristics based on the first time. ----------------------#
   model$npatches = npatches
   model$props    = data.frame( age   = pmax(1/4,ctrl$patch$age[[1]])
                              , lu    = ctrl$patch$lu           [[1]]
                              , lorey = ctrl$patch$can.depth    [[1]]
                              , agb   = ctrl$patch$agb          [[1]]
                              , lai   = ctrl$patch$lai          [[1]]
                              , ba    = ctrl$patch$ba           [[1]]
                              )#end list
   model$ptable   = ptable
   #---------------------------------------------------------------------------------------#


   #----- Find patches with the closest illumination to the quantiles to show. ------------#
   plist      = split(x=ptable,f=ptable$ipa)
   wfbeam          = mapply( FUN = function(x){
                                      ans = weighted.mean( x     = x$fbeam
                                                         , w     = x$cci.area/x$orig.area
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
                                      return(ans)
                                   }#end function
                           , x   = plist
                           )#end mapply
   wfarea          = mapply( FUN = function(x){
                                      ans = weighted.mean( x     = x$orig.area
                                                         , w     = x$cci.area/x$orig.area
                                                         , na.rm = TRUE
                                                         )#end weighted.mean
                                      return(ans)
                                   }#end function
                           , x   = plist
                           )#end mapply
   qfbeam          = quantile(x=wfbeam,probs=qshow,na.rm=TRUE)
   qidx.show       = mapply(FUN=which.closest,x=qfbeam,MoreArgs=list(A=wfbeam))
   model$qidx.show = qidx.show
   #---------------------------------------------------------------------------------------#



   #----- Find cohort, patch, and light indices so we store data as arrays. ---------------#
   ncohorts = 0
   for (e in eloop.ctrl){
      #----- Get time info. ---------------------------------------------------------------#
      now    = model$tomonth[e]
      mm     = nummonths(now)
      yyyy   = numyears (now)
      stamp  = paste0("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm))
      ncohorts = max(c(ncohorts,ctrl$cohort$ico[[stamp]]))
      #------------------------------------------------------------------------------------#
   }#end for (e in eloop.ctrl)
   for (e in eloop.test){
      #----- Get time info. ---------------------------------------------------------------#
      now    = model$tomonth[e]
      mm     = nummonths(now)
      yyyy   = numyears (now)
      stamp  = paste0("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm))
      ncohorts = max(c(ncohorts,test$cohort$ico[[stamp]]))
      #------------------------------------------------------------------------------------#
   }#end for (e in eloop.ctrl)
   model$ncohorts = ncohorts
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Define matrix placeholder.                                                        #
   #---------------------------------------------------------------------------------------#
   e.empty = array( data     = NA
                  , dim      = c(nemean,npatches,nlight)
                  , dimnames = list( as.character(model$tomonth)
                                   , sprintf("P%3.3i",sequence(npatches))
                                   , light.key
                                   )#end list
                  )#end array
   c.empty = array( data     = NA_real_
                  , dim      = c(nemean,3,2)
                  , dimnames = list( as.character(model$tomonth)
                                   , c("qlwr","mean","qupr")
                                   , c("ctrl","test")
                                   )#end list
                  )#end array
   q.empty = array( data     = NA_real_
                  , dim      = c(nemean,ndcycle,npatches,nlight)
                  , dimnames = list( as.character(model$tomonth)
                                   , dclabel
                                     , sprintf("P%3.3i",sequence(npatches))
                                   , light.key
                                   )#end list
                  )#end array
   l.empty = array( data     = NA_real_
                  , dim      = c(nemean,2,ncohorts,npatches,nlight)
                  , dimnames = list( as.character(model$tomonth)
                                   , c("cumlai","value")
                                   , sprintf("C%3.3i",sequence(ncohorts))
                                   , sprintf("P%3.3i",sequence(npatches))
                                   , light.key
                                   )#end list
                  )#end array
   hp.empty = array( data     = NA_real_
                   , dim      = c(nemean,n.hgt.lyr,3,2)
                   , dimnames = list( as.character(model$tomonth)
                                    , hgt.lyr.labels
                                    , c("qlwr","mean","qupr")
                                    , c("ctrl","test")
                                    )#end list
                   )#end array
   hc.empty = array( data     = NA_real_
                   , dim      = c(nemean,n.hgt.cdf,3,2)
                   , dimnames = list( as.character(model$tomonth)
                                    , hgt.cdf.labels
                                    , c("qlwr","mean","qupr")
                                    , c("ctrl","test")
                                    )#end list
                   )#end array
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Load all variables, and append them to matrices.                                 #
   #---------------------------------------------------------------------------------------#
   cat0("       ~ Load variables.")
   for (v in sequence(ncompvar)){
      #----- Load information. ------------------------------------------------------------#
      this.vnam      = compvar$vnam   [v]
      this.desc      = compvar$desc   [v]
      this.unit      = compvar$unit   [v]
      this.qmean     = compvar$qmean  [v]
      this.clprof    = compvar$clprof [v]
      this.clcum     = compvar$clcum  [v]
      this.hgtprof   = compvar$hgtprof[v]
      cat0("         > ",this.desc,".")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Initialise array for this variable.                                            #
      #------------------------------------------------------------------------------------#
      model[[this.vnam]] = list( emean = e.empty
                               , qmean = q.empty
                               , cmean = c.empty
                               , lprof = l.empty
                               , hprof = hp.empty
                               , cprof = hc.empty
                               )#end list
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Grab all formats of this variable.                                             #
      #------------------------------------------------------------------------------------#
      if (this.vnam %in% "tot.soil.c"){
         cpatch   = with( data = ctrl$patch
                        , expr = mapply( FUN      = sum.soil.c.fun
                                       , fast     = fast.soil.c
                                       , slow     = slow.soil.c
                                       , struct   = struct.soil.c
                                       , SIMPLIFY = FALSE
                                       )#end mapply
                        )#end with
         tpatch   = with( data = test$patch
                        , expr = mapply( FUN      = sum.soil.c.fun
                                       , fast     = fast.soil.c
                                       , slow     = slow.soil.c
                                       , struct   = struct.soil.c
                                       , SIMPLIFY = FALSE
                                       )#end mapply
                        )#end with
         q.cpatch = NULL
         q.tpatch = NULL
         #---------------------------------------------------------------------------------#
      }else if (this.vnam %in% "bowen"){
         cpatch   = with( data = ctrl$patch
                        , expr = mapply(FUN=bowen.fun,h=hflxca,qw=qwflxca,SIMPLIFY=FALSE)
                        )#end with
         tpatch   = with( data = test$patch
                        , expr = mapply(FUN=bowen.fun,h=hflxca,qw=qwflxca,SIMPLIFY=FALSE)
                        )#end with
         q.cpatch = with( data = ctrl$qpatch
                        , expr = mapply(FUN=bowen.fun,h=hflxca,qw=qwflxca,SIMPLIFY=FALSE)
                        )#end with
         q.tpatch = with( data = test$qpatch
                        , expr = mapply(FUN=bowen.fun,h=hflxca,qw=qwflxca,SIMPLIFY=FALSE)
                        )#end with
      }else if (this.vnam %in% "tratio"){
         cpatch   = with( data = ctrl$patch
                        , expr = mapply(FUN=tratio.fun,tp=transp,w=wflxca,SIMPLIFY=FALSE)
                       )#end with
         tpatch   = with( data = test$patch
                        , expr = mapply(FUN=tratio.fun,tp=transp,w=wflxca,SIMPLIFY=FALSE)
                       )#end with
         q.cpatch = with( data = ctrl$qpatch
                        , expr = mapply(FUN=tratio.fun,tp=transp,w=wflxca,SIMPLIFY=FALSE)
                        )#end with
         q.tpatch = with( data = test$qpatch
                        , expr = mapply(FUN=tratio.fun,tp=transp,w=wflxca,SIMPLIFY=FALSE)
                       )#end with
      }else if (this.vnam %in% c("zupr.gpp","zmid.gpp","zlwr.gpp")){
         args.top = modifyList(x=list(EXPR=substring(this.vnam,1,4)),val=clai.top.lim)
         args.bot = modifyList(x=list(EXPR=substring(this.vnam,1,4)),val=clai.bot.lim)
      
         ltop     = do.call(what="switch",args=args.top)
         lbot     = do.call(what="switch",args=args.bot)
         cpatch   = with( data     = ctrl$cohort
                        , expr     = mapply( FUN      = layer.gpp.list.fun
                                           , nplant   = nplant
                                           , area     = area
                                           , gpp      = gpp
                                           , lai      = lai
                                           , ipa      = ipa
                                           , MoreArgs = list(top=ltop,bot=lbot)
                                           , SIMPLIFY = FALSE
                                           )#end mapply
                        )#end with
         tpatch   = with( data     = test$cohort
                        , expr     = mapply( FUN      = layer.gpp.list.fun
                                           , nplant   = nplant
                                           , area     = area
                                           , gpp      = gpp
                                           , lai      = lai
                                           , ipa      = ipa
                                           , MoreArgs = list(top=ltop,bot=lbot)
                                           , SIMPLIFY = FALSE
                                           )#end mapply
                        )#end with
         q.cpatch = NULL
         q.tpatch = NULL
      }else{
         cpatch   = ctrl$patch [[this.vnam]]
         tpatch   = test$patch [[this.vnam]]
         q.cpatch = ctrl$qpatch[[this.vnam]]
         q.tpatch = test$qpatch[[this.vnam]]
      }#end if (this.vnam %in% "tot.soil.c")
      carea = ctrl$patch$area
      tarea = test$patch$area
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Load cohort-level variables.                                                   #
      #------------------------------------------------------------------------------------#
      if (is.na(this.clprof)){
         ccohort = NULL
         tcohort = NULL
      }else{
         ccohort = ctrl$cohort[[this.clprof]]
         tcohort = test$cohort[[this.clprof]]
      }#end if (is.na(this.clprof))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Copy control data to array.                                                    #
      #------------------------------------------------------------------------------------#
      for (e in eloop.ctrl){
         #----- Get time info. ------------------------------------------------------------#
         now    = model$tomonth[e]
         mm     = nummonths(now)
         yyyy   = numyears (now)
         stamp  = paste0("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm))
         #---------------------------------------------------------------------------------#



         #----- Copy homogeneous patch. ---------------------------------------------------#
         if(!is.null(cpatch  )) model[[this.vnam]]$emean[e, ,nlight] = cpatch[[stamp]]
         if(!is.null(q.cpatch)) model[[this.vnam]]$qmean[e,,,nlight] = t(q.cpatch[[stamp]])
         #---------------------------------------------------------------------------------#


         #----- Find mean and the lower and upper bounds of patch distribution. -----------#
         if (! is.null(cpatch)){
            #----- Find mean and the lower and the upper bounds. --------------------------#
            cmean = weighted.mean    (x=cpatch[[stamp]],w=carea[[stamp]]        )
            if (cband > 0){
               cqlwr = weighted.quantile(x=cpatch[[stamp]],w=carea[[stamp]],qu=qlwr)
               cqupr = weighted.quantile(x=cpatch[[stamp]],w=carea[[stamp]],qu=qupr)
               model[[this.vnam]]$cmean[e,,"ctrl"] = c(cqlwr,cmean,cqupr)
            }else{
               csdev = weighted.sd(x=cpatch[[stamp]],w=carea[[stamp]],M=npat.orig)
               model[[this.vnam]]$cmean[e,,"ctrl"] = cmean + c(-1.,0.,1.)*abs(cband)*csdev
            }#end if (cband > 0)
            #------------------------------------------------------------------------------#
         }#end if (! is.null(cpatch))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      In case cohort data exist, copy them to the lprof.                         #
         #---------------------------------------------------------------------------------#
         if (! is.null(ccohort)){
            #----- Find cohort, patch, and light indices so we store data as arrays. ------#
            c.ipa           = ctrl$cohort$ipa   [[stamp]]
            c.ilt           = nlight + 0L * c.ipa
            c.ico           = ctrl$cohort$ico   [[stamp]]
            c.idx           = cbind(c.ico,c.ipa,c.ilt)
            c.npa           = length(ctrl$patch$area[[stamp]])
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Cohort properties: keep both ctrl and test.                              #
            #------------------------------------------------------------------------------#
            c.pft           = ctrl$cohort$pft   [[stamp]]
            c.lai           = ctrl$cohort$lai   [[stamp]]
            c.hgt           = ctrl$cohort$height[[stamp]]
            #------------------------------------------------------------------------------#


            #------ Find cohort height layer and cdf quadrat. -----------------------------#
            c.ilyr = cut(x=c.hgt,breaks=hgt.lyr.brks,labels=hgt.lyr.labels,right=FALSE)
            c.icdf = cut(x=c.hgt,breaks=hgt.cdf.brks,labels=hgt.cdf.labels,right=FALSE)
            #------------------------------------------------------------------------------#



            #------ Find the cumulative LAI for each profile. -----------------------------#
            cidx                           = cbind(e,1,c.idx)
            model[[this.vnam]]$lprof[cidx] = c.lai
            #------------------------------------------------------------------------------#



            #----- Handy aliases. ---------------------------------------------------------#
            cvalueco  = ccohort[[stamp]]
            cidx      = cbind(e,2,c.idx)
            cresco    = c.lai > pft$lai.min[c.pft]
            if (is.na(this.clcum)){
               cscalco                        = rep(1.,times=length(cvalueco))
            }else if (this.clcum %in% "nplant"){
               cscalco = ctrl$cohort[[this.clcum]][[stamp]] / ctrl$cohort$area[[stamp]]
            }else{
               cscalco = ctrl$cohort[[this.clcum]][[stamp]]
            }#end if (is.na(this.clcum))
            #------------------------------------------------------------------------------#



            #----- Scale variable. --------------------------------------------------------#
            cintegco = ifelse(test=cresco,yes=cvalueco*cscalco,no=0.)
            careaco  = ctrl$cohort$area[[stamp]]
            #------------------------------------------------------------------------------#



            #----- Save the LAI profile. --------------------------------------------------#
            model[[this.vnam]]$lprof[cidx] = cintegco
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Find vertical profile at fixed height bins.                              #
            #------------------------------------------------------------------------------#
            if (this.hgtprof){

               #----- Find coarse-layer aggregated properties. ----------------------------#
               cinteg.agg      = tapply( X     = cintegco * careaco
                                       , INDEX = list(c.ilyr,c.ipa)
                                       , FUN   = sum
                                       )#end tapply
               cinteg.lyr      = matrix(data = 0,nrow=n.hgt.lyr,ncol=c.npa)
               carea.lyr       = matrix( data  = carea[[stamp]]
                                       , nrow  = n.hgt.lyr
                                       , ncol  = c.npa
                                       , byrow = TRUE
                                       )#end matrix
               idx             = cbind( match( rownames(cinteg.agg)[row(cinteg.agg)]
                                             , hgt.lyr.labels)
                                      , as.numeric(colnames(cinteg.agg)[col(cinteg.agg)])
                                      )#end cbind
               cinteg.lyr[idx] = ifelse(test=is.finite(cinteg.agg),yes=cinteg.agg,no=0.)
               cinteg.lyr      = split(x=cinteg.lyr,f=row(cinteg.lyr))
               carea.lyr       = split(x=carea.lyr,f=row(carea.lyr))
               cmean.lyr       = mapply( FUN      = weighted.mean
                                       , x        = cinteg.lyr
                                       , w        = carea.lyr
                                       )#end mapply
               if (cband > 0){
                  cqlwr.lyr    = mapply( FUN      = weighted.quantile
                                       , x        = cinteg.lyr
                                       , w        = carea.lyr
                                       , MoreArgs = list(qu=qlwr)
                                       )#end mapply
                  cqupr.lyr    = mapply( FUN      = weighted.quantile
                                       , x        = cinteg.lyr
                                       , w        = carea.lyr
                                       , MoreArgs = list(qu=qupr)
                                       )#end mapply
               }else{
                  csdev.lyr    = mapply( FUN      = weighted.sd
                                       , x        = cinteg.lyr
                                       , w        = carea.lyr
                                       , MoreArgs = list(M=npat.orig)
                                       )#end mapply
                  cqlwr.lyr    = cmean.lyr - 1. * abs(cband)*csdev.lyr
                  cqupr.lyr    = cmean.lyr + 1. * abs(cband)*csdev.lyr
               }#end if (cband > 0)
               cinteg.lyr = cbind(qlwr=cqlwr.lyr,mean=cmean.lyr,qupr=cqupr.lyr)
               #---------------------------------------------------------------------------#



               #----- Find fine-layer aggregated properties (with CDF). -------------------#
               cinteg.agg      = tapply(X=cintegco,INDEX=list(c.icdf,c.ipa),FUN=sum )
               cinteg.cdf      = matrix(data = 0,nrow=n.hgt.cdf,ncol=c.npa)
               carea.cdf       = matrix( data  = carea[[stamp]]
                                       , nrow  = n.hgt.cdf
                                       , ncol  = c.npa
                                       , byrow = TRUE
                                       )#end matrix
               idx             = cbind( match( rownames(cinteg.agg)[row(cinteg.agg)]
                                             , hgt.cdf.labels)
                                      , as.numeric(colnames(cinteg.agg)[col(cinteg.agg)])
                                      )#end cbind
               cinteg.cdf[idx] = ifelse(test=is.finite(cinteg.agg),yes=cinteg.agg,no=0.)
               cinteg.cdf      = apply( X      = cinteg.cdf
                                      , MARGIN = 2
                                      , FUN    = function(x) rev(cumsum(rev(x)))
                                      )#end apply
               cinteg.cdf      = split(x=cinteg.cdf,f=row(cinteg.cdf))
               carea.cdf       = split(x=carea.cdf ,f=row(carea.cdf))
               cmean.cdf       = mapply( FUN      = weighted.mean
                                       , x        = cinteg.cdf
                                       , w        = carea.cdf
                                       )#end mapply
               if (cband > 0){
                  cqlwr.cdf    = mapply( FUN      = weighted.quantile
                                       , x        = cinteg.cdf
                                       , w        = carea.cdf
                                       , MoreArgs = list(qu=qlwr)
                                       )#end mapply
                  cqupr.cdf    = mapply( FUN      = weighted.quantile
                                       , x        = cinteg.cdf
                                       , w        = carea.cdf
                                       , MoreArgs = list(qu=qupr)
                                       )#end mapply
               }else{
                  csdev.cdf    = mapply( FUN      = weighted.sd
                                       , x        = cinteg.cdf
                                       , w        = carea.cdf
                                       , MoreArgs = list(M=npat.orig)
                                       )#end mapply
                  cqlwr.cdf    = cmean.cdf - 1. * abs(cband)*csdev.cdf
                  cqupr.cdf    = cmean.cdf + 1. * abs(cband)*csdev.cdf
               }#end if (cband > 0)
               cinteg.cdf = cbind(qlwr=cqlwr.cdf,mean=cmean.cdf,qupr=cqupr.cdf)
               #---------------------------------------------------------------------------#


               #----- Copy results. -------------------------------------------------------#
               model[[this.vnam]]$hprof[e,,"qlwr","ctrl"] = cinteg.lyr[,"qlwr"]
               model[[this.vnam]]$hprof[e,,"mean","ctrl"] = cinteg.lyr[,"mean"]
               model[[this.vnam]]$hprof[e,,"qupr","ctrl"] = cinteg.lyr[,"qupr"]
               model[[this.vnam]]$cprof[e,,"qlwr","ctrl"] = cinteg.cdf[,"qlwr"]
               model[[this.vnam]]$cprof[e,,"mean","ctrl"] = cinteg.cdf[,"mean"]
               model[[this.vnam]]$cprof[e,,"qupr","ctrl"] = cinteg.cdf[,"qupr"]
               #---------------------------------------------------------------------------#
            }#end if (this.hgtprof)
            #------------------------------------------------------------------------------#
         }#end if (! is.null(ccohort))
         #---------------------------------------------------------------------------------#
      }#end for (e in emean.loop)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Copy test data to array.                                                       #
      #------------------------------------------------------------------------------------#
      for (e in eloop.test){
         #----- Get time info. ------------------------------------------------------------#
         now   = model$tomonth[e]
         mm    = nummonths(now)
         yyyy  = numyears (now)
         stamp = paste0("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Make sure variable exists.                                                  #
         #---------------------------------------------------------------------------------#
         if (! is.null(tpatch)){
            tnow  = tpatch[[stamp]]

            #----- Copy data for each light type (some may be missing). -------------------#
            model[[this.vnam]]$emean[e,,nwl] = 0.
            for (l in sequence(nltype)){
               sel = ptable$ilight == l
               idx = ptable$ipa[sel]
               wgt = ptable$cci.area[sel] / ptable$orig.area[sel]
               model[[this.vnam]]$emean[e,idx,l]   = tnow[sel]
               model[[this.vnam]]$emean[e,idx,nwl] = ( model[[this.vnam]]$emean[e,idx,nwl]
                                                     + tnow[sel] * wgt
                                                     )#end
            }#end for (l in sequence(nltype))
            #------------------------------------------------------------------------------#


            #----- Find mean and the lower and the upper bounds. --------------------------#
            tmean = weighted.mean    (x=tpatch[[stamp]],w=tarea[[stamp]]        )
            if (cband > 0){
               tqlwr = weighted.quantile(x=tpatch[[stamp]],w=tarea[[stamp]],qu=qlwr)
               tqupr = weighted.quantile(x=tpatch[[stamp]],w=tarea[[stamp]],qu=qupr)
               model[[this.vnam]]$cmean[e,,"test"] = c(tqlwr,tmean,tqupr)
            }else{
               tsdev = weighted.sd(x=tpatch[[stamp]],w=tarea[[stamp]],M=npat.orig)
               model[[this.vnam]]$cmean[e,,"test"] = tmean + c(-1.,0.,1.)*abs(cband)*tsdev
            }#end if (cband > 0)
            #------------------------------------------------------------------------------#
         }#end if (! is.null(tpatch))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Make sure variable exists.                                                  #
         #---------------------------------------------------------------------------------#
         if (! is.null(q.tpatch)){
            qnow  = t(q.tpatch[[stamp]])
            #----- Copy data for each light type (some may be missing). -------------------#
            model[[this.vnam]]$qmean[e,,,nwl] = 0.
            for (l in sequence(nltype)){
               sel = ptable$ilight == l
               idx = ptable$ipa[sel]
               wgt = ptable$cci.area[sel] / ptable$orig.area[sel]
               wgt = matrix(wgt,nrow=ndcycle,ncol=length(wgt),byrow=TRUE)
               model[[this.vnam]]$qmean[e,,idx,l ]  = qnow[,sel]
               model[[this.vnam]]$qmean[e,,idx,nwl] = ( model[[this.vnam]]$qmean[e,,idx,nwl]
                                                      + qnow[,sel] * wgt
                                                      )#end
            }#end for (l in sequence(nltype))
            #------------------------------------------------------------------------------#
         }#end if (! is.null(tpatch))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      In case cohort data exist, copy them to the lprof.                         #
         #---------------------------------------------------------------------------------#
         if (! is.null(tcohort)){
            #----- Find cohort, patch, and light indices so we store data as arrays. ------#
            t.ipa           = ptable$ipa   [test$cohort$ipa[[stamp]]]
            t.ilt           = ptable$ilight[test$cohort$ipa[[stamp]]]
            t.ico           = test$cohort$ico   [[stamp]]
            t.idx           = cbind(t.ico,t.ipa,t.ilt)
            t.npa           = length(test$patch$area[[stamp]])
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Cohort properties: keep both test and test.                              #
            #------------------------------------------------------------------------------#
            t.pft           = test$cohort$pft   [[stamp]]
            t.lai           = test$cohort$lai   [[stamp]]
            t.hgt           = test$cohort$height[[stamp]]
            #------------------------------------------------------------------------------#


            #------ Find cohort height layer and cdf quadrat. -----------------------------#
            t.ilyr = cut(x=t.hgt,breaks=hgt.lyr.brks,labels=hgt.lyr.labels,right=FALSE)
            t.icdf = cut(x=t.hgt,breaks=hgt.cdf.brks,labels=hgt.cdf.labels,right=FALSE)
            #------------------------------------------------------------------------------#



            #------ Find the cumulative LAI for each profile. -----------------------------#
            tidx                           = cbind(e,1,t.idx)
            model[[this.vnam]]$lprof[tidx] = t.lai
            #------------------------------------------------------------------------------#


            #----- Handy aliases. ---------------------------------------------------------#
            tvalueco  = tcohort[[stamp]]
            tidx      = cbind(e,2,t.idx)
            tresco    = t.lai > pft$lai.min[t.pft]
            if (is.na(this.clcum)){
               tscalco                        = rep(1.,times=length(tvalueco))
            }else if (this.clcum %in% "nplant"){
               tscalco = test$cohort[[this.clcum]][[stamp]] / test$cohort$area[[stamp]]
            }else{
               tscalco = test$cohort[[this.clcum]][[stamp]]
            }#end if (is.na(this.clcum))
            #------------------------------------------------------------------------------#



            #----- Scale variable. --------------------------------------------------------#
            tintegco = ifelse(test=tresco,yes=tvalueco*tscalco,no=0.)
            tareaco  = test$cohort$area[[stamp]]
            #------------------------------------------------------------------------------#


            #----- Save the LAI profile. --------------------------------------------------#
            model[[this.vnam]]$lprof[tidx] = tintegco
            #------------------------------------------------------------------------------#



            #----- Find weighted average. -------------------------------------------------#
            wgt      = matrix(data=0,ncol=nltype,nrow=npatches)
            idx      = cbind(ptable$ipa,ptable$ilight)
            wgt[idx] = ptable$cci.area / ptable$orig.area
            lprof    = model[[this.vnam]]$lprof
            for (ip in sequence(npatches)){
               for (ix in sequence(2)){
                  lprof[e,ix,,ip,nwl] = apply( X      = lprof[e,ix,,ip,sequence(nltype)]
                                             , MARGIN = 1
                                             , FUN    = weighted.mean
                                             , w      = wgt[ip,]
                                             , na.rm  = TRUE
                                             )#end apply
               }#end for (ix in sequence(2))
            }#end for (ip in sequence(npatches))
            model[[this.vnam]]$lprof = ifelse(is.finite(lprof),lprof,NA)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find vertical profile at fixed height bins.                              #
            #------------------------------------------------------------------------------#
            if (this.hgtprof){

               #----- Find coarse-layer aggregated properties. ----------------------------#
               tinteg.agg      = tapply( X     = tintegco*tareaco
                                       , INDEX = list(t.ilyr,t.ipa)
                                       , FUN   = sum
                                       )#end tapply
               tinteg.lyr      = matrix(data = 0,nrow=n.hgt.lyr,ncol=t.npa)
               tarea.lyr       = matrix( data  = tarea[[stamp]]
                                       , nrow  = n.hgt.lyr
                                       , ncol  = t.npa
                                       , byrow = TRUE
                                       )#end matrix
               idx             = cbind( match( rownames(tinteg.agg)[row(tinteg.agg)]
                                             , hgt.lyr.labels)
                                      , as.numeric(colnames(tinteg.agg)[col(tinteg.agg)])
                                      )#end cbind
               tinteg.lyr[idx] = ifelse(test=is.finite(tinteg.agg),yes=tinteg.agg,no=0.)
               tinteg.lyr      = split(x=tinteg.lyr,f=row(tinteg.lyr))
               tarea.lyr       = split(x=tarea.lyr,f=row(tarea.lyr))
               tmean.lyr       = mapply( FUN      = weighted.mean
                                       , x        = tinteg.lyr
                                       , w        = tarea.lyr
                                       )#end mapply
               if (cband > 0){
                  tqlwr.lyr    = mapply( FUN      = weighted.quantile
                                       , x        = tinteg.lyr
                                       , w        = tarea.lyr
                                       , MoreArgs = list(qu=qlwr)
                                       )#end mapply
                  tqupr.lyr    = mapply( FUN      = weighted.quantile
                                       , x        = tinteg.lyr
                                       , w        = tarea.lyr
                                       , MoreArgs = list(qu=qupr)
                                       )#end mapply
               }else{
                  tsdev.lyr    = mapply( FUN      = weighted.sd
                                       , x        = tinteg.lyr
                                       , w        = tarea.lyr
                                       , MoreArgs = list(M=npat.orig)
                                       )#end mapply
                  tqlwr.lyr    = tmean.lyr - 1. * abs(cband)*tsdev.lyr
                  tqupr.lyr    = tmean.lyr + 1. * abs(cband)*tsdev.lyr
               }#end if (cband > 0)
               tinteg.lyr = cbind(qlwr=tqlwr.lyr,mean=tmean.lyr,qupr=tqupr.lyr)
               #---------------------------------------------------------------------------#



               #----- Find fine-layer aggregated properties (with CDF). -------------------#
               tinteg.agg      = tapply(X=tintegco,INDEX=list(t.icdf,t.ipa),FUN=sum )
               tinteg.cdf      = matrix(data = 0,nrow=n.hgt.cdf,ncol=t.npa)
               tarea.cdf       = matrix( data  = tarea[[stamp]]
                                       , nrow  = n.hgt.cdf
                                       , ncol  = t.npa
                                       , byrow = TRUE
                                       )#end matrix
               idx             = cbind( match( rownames(tinteg.agg)[row(tinteg.agg)]
                                             , hgt.cdf.labels)
                                      , as.numeric(colnames(tinteg.agg)[col(tinteg.agg)])
                                      )#end cbind
               tinteg.cdf[idx] = ifelse(test=is.finite(tinteg.agg),yes=tinteg.agg,no=0.)
               tinteg.cdf      = apply( X      = tinteg.cdf
                                      , MARGIN = 2
                                      , FUN    = function(x) rev(cumsum(rev(x)))
                                      )#end apply
               tinteg.cdf      = split(x=tinteg.cdf,f=row(tinteg.cdf))
               tarea.cdf       = split(x=tarea.cdf ,f=row(tarea.cdf))
               tmean.cdf       = mapply( FUN      = weighted.mean
                                       , x        = tinteg.cdf
                                       , w        = tarea.cdf
                                       )#end mapply
               if (cband > 0){
                  tqlwr.cdf    = mapply( FUN      = weighted.quantile
                                       , x        = tinteg.cdf
                                       , w        = tarea.cdf
                                       , MoreArgs = list(qu=qlwr)
                                       )#end mapply
                  tqupr.cdf    = mapply( FUN      = weighted.quantile
                                       , x        = tinteg.cdf
                                       , w        = tarea.cdf
                                       , MoreArgs = list(qu=qupr)
                                       )#end mapply
               }else{
                  tsdev.cdf    = mapply( FUN      = weighted.sd
                                       , x        = tinteg.cdf
                                       , w        = tarea.cdf
                                       , MoreArgs = list(M=npat.orig)
                                       )#end mapply
                  tqlwr.cdf    = tmean.cdf - 1. * abs(cband)*tsdev.cdf
                  tqupr.cdf    = tmean.cdf + 1. * abs(cband)*tsdev.cdf
               }#end if (cband > 0)
               tinteg.cdf = cbind(qlwr=tqlwr.cdf,mean=tmean.cdf,qupr=tqupr.cdf)
               #---------------------------------------------------------------------------#


               #----- Copy results. -------------------------------------------------------#
               model[[this.vnam]]$hprof[e,,"qlwr","test"] = tinteg.lyr[,"qlwr"]
               model[[this.vnam]]$hprof[e,,"mean","test"] = tinteg.lyr[,"mean"]
               model[[this.vnam]]$hprof[e,,"qupr","test"] = tinteg.lyr[,"qupr"]
               model[[this.vnam]]$cprof[e,,"qlwr","test"] = tinteg.cdf[,"qlwr"]
               model[[this.vnam]]$cprof[e,,"mean","test"] = tinteg.cdf[,"mean"]
               model[[this.vnam]]$cprof[e,,"qupr","test"] = tinteg.cdf[,"qupr"]
               #---------------------------------------------------------------------------#
            }#end if (this.hgtprof)
            #------------------------------------------------------------------------------#
         }#end if (! is.null(tcohort))
         #---------------------------------------------------------------------------------#
      }#end for (e in emean.loop)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Find the weighted mean for each layer of vertical profile.                    #
      #------------------------------------------------------------------------------------#
      if (! ( is.na(this.clprof) || is.na(this.clcum))){
         lprof                    = apply( X = model[[this.vnam]]$lprof
                                         , MARGIN = c(1,2,4,5)
                                         , FUN    = cumsum
                                         )#end apply
         model[[this.vnam]]$lprof = aperm(a=lprof,perm=c(2,3,1,4,5))
      }#end if (! ( is.na(this.clprof) || is.na(this.clcum)))
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#



   #----- Copy the data to the results. ---------------------------------------------------#
   res.iata    = paste("res",iata,sep=".")
   assign(res.iata,model)
   res[[iata]] = model
   rm(list=c("model"))
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      Save processed data to RData.                                                    #
   #---------------------------------------------------------------------------------------#
   rdata.iata = file.path(rdata.path,paste(iata,sim.suffix,rdata.suffix,sep="_"))
   cat0(" + Save processed data to ",basename(rdata.iata),".")
   dummy = save(list=c(res.iata), file=rdata.iata)
   rm(res.iata,this.sim)
   #---------------------------------------------------------------------------------------#

}#end for (s in sequence(loop.nsites))
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Plot the time series by patch.                                                       #
#------------------------------------------------------------------------------------------#
cat0(" + Plot long-term time series.")
for (s in sequence(nsites)){
   #----- Loop over variables. ------------------------------------------------------------#
   iata     = sites$iata[s]
   longname = sites$desc[s]
   model    = res[[iata]]
   npatches = model$npatches
   props    = model$props
   ptable   = model$ptable
   #---------------------------------------------------------------------------------------#
   
   cat0("   - ",longname,".")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.vnam       = compvar$vnam[v]
      this.desc       = compvar$desc[v]
      this.unit       = untab[[compvar$unit[v]]]
      this.qmean      = compvar$qmean[v]
      this.cscheme    = get(compvar$cscheme[v])
      this.clprof     = ! is.na(compvar$clprof[v])
      this.clcum      = ! is.na(compvar$clcum [v])
      this.hgtprof    = compvar$hgtprof[v]
      tomonth         = model$tomonth
      toyear          = sort(unique(model$toyear))
      emean           = model[[this.vnam]]$emean
      cmean           = model[[this.vnam]]$cmean
      ymean           = qapply(X=emean,INDEX=model$toyear,DIM=1,FUN=mean,na.rm=TRUE)
      qmean           = model[[this.vnam]]$qmean
      lprof           = model[[this.vnam]]$lprof
      wshow           = model$tomonth %wr% when.show
      qidx.show       = model$qidx.show
      cat0("     > ",this.desc,".")
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Decide whether to compute or create a dummy qmean array.
      #------------------------------------------------------------------------------------#
      if (this.qmean){
         qmean        = qapply( X     = qmean[wshow,,,,drop=FALSE]
                              , INDEX = nummonths(model$tomonth[wshow])
                              , DIM   = 1
                              , FUN   = mean
                              , na.rm = TRUE
                              )#end qapply
      }else{
         qmean        = array(data=NA,dim=c(12,model$ndcycle,dim(emean)[c(2,3)]))
      }#end if (this.qmean)
      dimnames(qmean) = list( month.abb
                            , model$dclabel
                            , dimnames(emean)[[2]]
                            , dimnames(emean)[[3]]
                            )#end list
      #------------------------------------------------------------------------------------#




      #----- Limit times for the profiles. ------------------------------------------------#
      hprof           = model[[this.vnam]]$hprof[wshow,,,,drop=FALSE]
      cprof           = model[[this.vnam]]$cprof[wshow,,,,drop=FALSE]
      hplab           = list(month.abb,dimnames(hprof)[[2]]
                            ,dimnames(hprof)[[3]],dimnames(hprof)[[4]])
      cplab           = list(month.abb,dimnames(cprof)[[2]]
                            ,dimnames(cprof)[[3]],dimnames(hprof)[[4]])
      hprof           = qapply( X     = hprof
                              , INDEX = nummonths(model$tomonth[wshow])
                              , DIM   = 1
                              , FUN   = mean
                              , na.rm = TRUE
                              )#end qapply
      cprof           = qapply( X     = cprof
                              , INDEX = nummonths(model$tomonth[wshow])
                              , DIM   = 1
                              , FUN   = mean
                              , na.rm = TRUE
                              )#end qapply
      dimnames(hprof) = hplab
      dimnames(cprof) = cplab
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #      Reduce the number of dimensions for vertical profile.                         #
      #------------------------------------------------------------------------------------#
      if (this.clprof){
         lprof        = qapply( X     = lprof[wshow,,,,,drop=FALSE]
                              , INDEX = nummonths(model$tomonth[wshow])
                              , DIM   = 1
                              , FUN   = mean
                              , na.rm = TRUE
                              )#end qapply
      }else{
         lprof        = array(data=NA,dim=c(12,2,model$ncohorts,dim(emean)[c(2,3)]))
      }#end if (this.qmean)
      dimnames(lprof) = list( month.abb
                            , c("cumlai","value")
                            , sprintf("C%3.3i",sequence(model$ncohorts))
                            , dimnames(emean)[[2]]
                            , dimnames(emean)[[3]]
                            )#end list
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find plot limits for time.                                                     #
      #------------------------------------------------------------------------------------#      
      em.xlimit = chron(range(x=model$tomonth[wshow]))
      sm.xlimit = em.xlimit
      ym.xlimit = pretty.xylim(u=model$toyear [wshow])
      qm.xlimit = pretty.xylim(u=model$dchour)
      em.pretty = pretty.time (when=em.xlimit)
      sm.pretty = em.pretty
      qm.pretty = pretty.elapsed(x=model$dchour,base=24)
      lm.yrange = apply(X=lprof[,"cumlai",,,,drop=FALSE],MARGIN=4,FUN=range,finite=TRUE)
      #------------------------------------------------------------------------------------#

      #------ Limits for height profiles. -------------------------------------------------#
      hm.xlimit = pretty.xylim(u=c(hprof))
      hm.xpretty = pretty(hm.xlimit)
      hm.xlabels = sprintf("%g",hm.xpretty)
      cm.xlimit = pretty.xylim(u=c(cprof))
      cm.xpretty = pretty(cm.xlimit)
      cm.xlabels = sprintf("%g",cm.xpretty)
      hm.ylimit  = pretty.xylim(u=hgt.lyr.range)
      hm.ypretty = pretty(hm.ylimit)
      hm.ylabels = sprintf("%g",hm.ypretty)
      cm.ylimit  = pretty.xylim(u=hgt.cdf.range)
      cm.ypretty = pretty(cm.ylimit)
      cm.ylabels = sprintf("%g",cm.ypretty)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find range for each patch.                                                     #
      #------------------------------------------------------------------------------------#
      eshow     = chron(dimnames(emean)[[1]])      %wr% when.show
      cshow     = chron(dimnames(cmean)[[1]])      %wr% when.show
      yshow     = as.numeric(dimnames(ymean)[[1]]) %wr% numyears(when.show)
      em.yrange = apply(X=emean[eshow,,,drop=FALSE]    ,MARGIN=2,FUN=range,finite=TRUE)
      sm.yrange = range(c(cmean[eshow,,]),finite=TRUE)
      ym.yrange = apply(X=ymean[yshow,,,drop=FALSE]    ,MARGIN=2,FUN=range,finite=TRUE)
      qm.yrange = apply(X=qmean                        ,MARGIN=3,FUN=range,finite=TRUE)
      mm.yrange = apply(X=emean                        ,MARGIN=1,FUN=range,finite=TRUE)
      lm.xrange = apply(X=lprof[,"value",,,,drop=FALSE],MARGIN=4,FUN=range,finite=TRUE)
      #------------------------------------------------------------------------------------#




      #------ Set some common features. ---------------------------------------------------#
      leclai  = desc.unit(desc="Cumulative LAI",unit=untab$m2lom2)
      ley     = desc.unit(desc=this.desc       ,unit=this.unit   )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over patches.                                                             #
      #------------------------------------------------------------------------------------#
      for (iq in sequence(nquant)){
         #----- Assign label and suffix. --------------------------------------------------#
         p         = qidx.show [iq]
         qkey      = qkey.show [iq]
         qdesc     = qdesc.show[iq]
         qsuffix   = paste0(this.vnam,"-",iata,"_",qkey)
         outepatch = file.path(outemean,qkey    )
         outypatch = file.path(outymean,qkey    )
         if (! file.exists(outepatch)             ) dir.create(outepatch)
         if (! file.exists(outypatch) & plot.ymean) dir.create(outypatch)
         if (this.qmean){
            outqpmain = file.path(outqmean ,this.vnam)
            outqpatch = file.path(outqpmain,qkey     )
            if (! file.exists(outqpmain)) dir.create(outqpmain)
            if (! file.exists(outqpatch)) dir.create(outqpatch)
         }#end if (this.qmean)
         if (this.clprof){
            outclpmain = file.path(outcumlai,this.vnam)
            outclpatch = file.path(outclpmain,qkey    )
            if (! file.exists(outclpmain)) dir.create(outclpmain)
            if (! file.exists(outclpatch)) dir.create(outclpatch)
         }#end if (this.clprof)
         outhpmain = file.path(outhgtprof,this.vnam)
         if (! file.exists(outhpmain)) dir.create(outhpmain)
         outhcmain = file.path(outhgtcum,this.vnam)
         if (! file.exists(outhcmain)) dir.create(outhcmain)
         outmpmain = file.path(outmmean,this.vnam)
         if (! file.exists(outmpmain)) dir.create(outmpmain)
         #---------------------------------------------------------------------------------#


         #----- Limits for this patch. ----------------------------------------------------#
         em.ylimit  = pretty.xylim(u= c(em.yrange[,p]))
         ym.ylimit  = pretty.xylim(u= c(ym.yrange[,p]))
         qm.ylimit  = pretty.xylim(u= c(qm.yrange[,p]))
         lm.xlimit  = pretty.xylim(u= c(lm.xrange[,p]))
         lm.xpretty = pretty(lm.xlimit)
         lm.xlabels = sprintf("%g",lm.xpretty)
         lm.ylimit  = pretty.xylim(u=-c(lm.yrange[,p]))
         lm.ypretty = pretty(lm.ylimit)
         lm.ylabels = sprintf("%g",ifelse(lm.ypretty == 0.,0,-lm.ypretty))
         #---------------------------------------------------------------------------------#


         #----- Make title. ---------------------------------------------------------------#
         le.emean = paste0(longname,"\n","Monthly means - "  ,qdesc )
         le.ymean = paste0(longname,"\n","Annual means - "   ,qdesc )
         #---------------------------------------------------------------------------------#


         #----- Paste patch properties. ---------------------------------------------------#
         leage = paste0("Age = "           ,sprintf("%7.2f",props$age  [p]))
         lelai = paste0("LAI = "           ,sprintf("%7.2f",props$lai  [p]))
         lelor = paste0("Lorey's height = ",sprintf("%7.2f",props$lorey[p]))
         lelab = c( age = desc.unit(desc=leage,unit=untab$yr    ,bracket=FALSE)
                  , lai = desc.unit(desc=lelai,unit=untab$m2lom2,bracket=FALSE)
                  , lor = desc.unit(desc=lelor,unit=untab$m     ,bracket=FALSE)
                  )#end c
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Create detailed legend information.                                         #
         #---------------------------------------------------------------------------------#
         light.lab = rep(NA_character_,times=nlight)
         for (l in sequence(nltype)){
            sel = ptable$ilight == l & ptable$ipa == p
            if (any(sel)){
               a.frac       = paste0( "Area = "
                                    , sprintf("%.1f" ,100. * ptable$cci.area [sel]
                                                           / ptable$orig.area[sel]
                                             )#end sprintf
                                    , "%"
                                    )#end paste0
               a.light      = paste0( "Light = "
                                    , sprintf("%.1f",100.*ptable$fbeam[sel])
                                    , "%"
                                    )#end paste0
               light.lab[l] = paste0(light.key[l],"; ",a.frac,"; ",a.light)
            }else{
               light.lab[l] = paste0(light.key[l],"; Non-existent.")
            }#end if
         }#end for  (l in sequence(nlight-2))
         #---------------------------------------------------------------------------------#

         #---- Average illumination for this patch. ---------------------------------------#
         l            = nwl
         sel          = ptable$ipa == p
         a.light      = weighted.mean( x = ptable$fbeam[sel]
                                     , w = ptable$cci.area[sel] / ptable$orig.area[sel]
                                     )#end weighted.mean
         a.light      = paste0(" Avg. light = ",sprintf("%.1f",100*a.light),"%")
         light.lab[l] = paste0(light.key[l],"; ",a.light)
         #---------------------------------------------------------------------------------#



         #----- Fraction of landscape covered by the area. --------------------------------#
         l            = nlight
         sel          = ptable$ipa == p
         a.frac       = mean(x = ptable$orig.area[sel])
         a.light      = paste0(" Fraction of landscape = ",sprintf("%.1f",100*a.frac),"%")
         light.lab[l] = paste0(light.key[l],"; ",a.light)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Plot monthly means.                                                        #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            fichier = file.path(outepatch,paste0("emean-",qsuffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = xsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            layout(mat= rbind(c(3,3),c(1,2)),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #----- Plot patch properties. -------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = lelab
                  , fill    = "transparent"
                  , border  = "transparent"
                  , ncol    = 1
                  , cex     = 0.8
                  , xpd     = TRUE
                  , bty     = "n"
                  )#end legend
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = light.lab
                  , col     = light.col
                  , lty     = light.lty
                  , lwd     = light.lwd
                  , ncol    = 1
                  , cex     = 0.8
                  , xpd     = TRUE
                  , bty     = "n"
                  )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot monthly means.                                                     #
            #------------------------------------------------------------------------------#
            par(mar=c(3.1,4.6,3.1,1.6))
            plot.new()
            plot.window(xlim=em.xlimit,ylim=em.ylimit)
            axis(side=1,las=1,at=em.pretty$levels,labels=em.pretty$labels)
            axis(side=2,las=1)
            title(main=le.emean,ylab=ley,cex.main=1.0)
            for (l in sequence(nlight)){
               #----- Load variables. -----------------------------------------------------#
               lines( x    = tomonth[eshow]
                    , y    = emean  [eshow,p,l]
                    , col  = light.col[l]
                    , lwd  = light.lwd[l]
                    , lty  = light.lty[l]
                    , type = "l"
                    )#end lines
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            box()
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#






         #---------------------------------------------------------------------------------#
         #      Plot annual means.                                                         #
         #---------------------------------------------------------------------------------#
         if (plot.ymean){ymean.loop = sequence(nout)}else{ymean.loop = numeric(0)}
         for (o in ymean.loop){
            #----- Make the file name. ----------------------------------------------------#
            fichier = file.path(outypatch,paste0("ymean-",qsuffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = xsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            layout(mat= rbind(c(3,3),c(1,2)),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #----- Plot patch properties. -------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = lelab
                  , fill    = "transparent"
                  , border  = "transparent"
                  , ncol    = 1
                  , cex     = 0.8
                  , xpd     = TRUE
                  , bty     = "o"
                  )#end legend
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = light.lab
                  , col     = light.col
                  , lty     = light.lty
                  , lwd     = light.lwd
                  , ncol    = 1
                  , cex     = 0.8
                  , xpd     = TRUE
                  , bty     = "o"
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot annual means.                                                      #
            #------------------------------------------------------------------------------#
            par(mar=c(3.1,4.6,3.1,1.6))
            plot.new()
            plot.window(xlim=ym.xlimit,ylim=ym.ylimit)
            axis(side=1,las=1,at=toyear,labels=sprintf("%g",toyear))
            axis(side=2,las=1)
            title(main=le.ymean,ylab=ley,cex.main=1.0)
            for (l in sequence(nlight)){
               lines( x    = toyear[yshow]
                    , y    = ymean[yshow,p,l]
                    , col  = light.col[l]
                    , lwd  = light.lwd[l]
                    , lty  = light.lty[l]
                    , type = "l"
                    )#end lines
            }#end for (s in sequence(nsimul))
            box()
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#






         #---------------------------------------------------------------------------------#
         #      In case the mean diel is to be plotted for this variable, go through each  #
         # month.                                                                          #
         #---------------------------------------------------------------------------------#
         mon.loop = sequence(12*this.qmean)
         for (m in mon.loop){
            #----- Aliases for current month. ---------------------------------------------#
            mlabel   = paste0("M",sprintf("%2.2i",m))
            mdesc    = month.name[m]
            mqsuffix = paste0(this.vnam,"-",iata,"_",mlabel,"_",qkey)
            le.qmean = paste0(longname,"\n","Mean diurnal cycle (",mdesc,") ",qdesc)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot monthly mean diurnal cycle.                                        #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = file.path(outqpatch,paste0("qmean-",mqsuffix,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = xsize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               layout(mat= rbind(c(3,3),c(1,2)),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #----- Plot patch properties. ----------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = lelab
                     , fill    = "transparent"
                     , border  = "transparent"
                     , ncol    = 1
                     , cex     = 0.8
                     , xpd     = TRUE
                     , bty     = "o"
                     )#end legend
               #---------------------------------------------------------------------------#



               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = light.lab
                     , col     = light.col
                     , lty     = light.lty
                     , lwd     = light.lwd
                     , ncol    = 1
                     , cex     = 0.8
                     , xpd     = TRUE
                     , bty     = "o"
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot annual means.                                                   #
               #---------------------------------------------------------------------------#
               par(mar=c(3.1,4.6,3.1,1.6))
               plot.new()
               plot.window(xlim=qm.xlimit,ylim=qm.ylimit)
               axis(side=1,las=1,at=qm.pretty,labels=sprintf("%2.2i",qm.pretty))
               axis(side=2,las=1)
               title(main=le.qmean,ylab=ley,cex.main=1.0)
               for (l in sequence(nlight)){
                  lines( x    = model$dchour
                       , y    = qmean[m,,p,l]
                       , col  = light.col[l]
                       , lwd  = light.lwd[l]
                       , lty  = light.lty[l]
                       , type = "l"
                       )#end lines
               }#end for (s in sequence(nsimul))
               box()
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end for (m in mon.loop)
         #---------------------------------------------------------------------------------#






         #---------------------------------------------------------------------------------#
         #      In case the mean profile is to be plotted for this variable, go through    #
         # each month.                                                                     #
         #---------------------------------------------------------------------------------#
         mon.loop = sequence(12*this.clprof)
         for (m in mon.loop){
            #----- Aliases for current month. ---------------------------------------------#
            mlabel   = paste0("M",sprintf("%2.2i",m))
            mdesc    = month.name[m]
            mzsuffix = paste0(this.vnam,"-",iata,"_",mlabel,"_",qkey)
            if (this.clcum){
               le.lprof = paste0(longname,"\n","Cumulative profile (",mdesc,") ",qdesc)
            }else{
               le.lprof = paste0(longname,"\n","Vertical profile (",mdesc,") ",qdesc)
            }#end if (this.clcum)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot monthly mean diurnal cycle.                                        #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = file.path(outclpatch,paste0("lprof-",mzsuffix,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = xsize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               layout(mat= rbind(c(3,3),c(1,2)),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #----- Plot patch properties. ----------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = lelab
                     , fill    = "transparent"
                     , border  = "transparent"
                     , ncol    = 1
                     , cex     = 0.8
                     , xpd     = TRUE
                     , bty     = "o"
                     )#end legend
               #---------------------------------------------------------------------------#



               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = light.lab
                     , col     = light.col
                     , lty     = light.lty
                     , lwd     = light.lwd
                     , ncol    = 1
                     , cex     = 0.8
                     , xpd     = TRUE
                     , bty     = "o"
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot annual means.                                                   #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.6,3.1,1.6))
               plot.new()
               plot.window(xlim=lm.xlimit,ylim=lm.ylimit)
               axis(side=1,las=1,at=lm.xpretty,labels=lm.xlabels)
               axis(side=2,las=1,at=lm.ypretty,labels=lm.ylabels)
               title(main=le.lprof,xlab=ley,ylab=leclai,cex.main=1.0)
               for (l in sequence(nlight)){
                  sel = is.finite(lprof[m,2,,p,l])
                  lines( x    = lprof[m,2,sel,p,l]
                       , y    = -lprof[m,1,sel,p,l]
                       , col  = light.col[l]
                       , lwd  = light.lwd[l]
                       , lty  = light.lty[l]
                       , type = "l"
                       )#end lines
               }#end for (s in sequence(nsimul))
               box()
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end for (m in mon.loop)
         #---------------------------------------------------------------------------------#

      }#end for (p in sequence(npatches))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot mean and horizontal heterogeneity.                                        #
      #------------------------------------------------------------------------------------#
         #----- Assign label and suffix. --------------------------------------------------#
         csuffix   = paste0(this.vnam,"-",iata)
         #---------------------------------------------------------------------------------#



         #----- Limits for this patch. ----------------------------------------------------#
         sm.ylimit = pretty.xylim(u=c(sm.yrange    ))
         #---------------------------------------------------------------------------------#


         #----- Make title. ---------------------------------------------------------------#
         le.cmean = paste0(longname,"\n","Monthly means and ",cbdesc)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Plot monthly means.                                                        #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            fichier = file.path(outcmean,paste0("cmean-",csuffix,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = xsize
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            layout(mat= rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = cband.leg
                  , fill    = cband.col
                  , border  = cband.col
                  , angle   = cband.angle
                  , density = cband.density
                  , col     = cmean.col
                  , lty     = cmean.lty
                  , lwd     = 2.5
                  , ncol    = 2
                  , cex     = 1.0
                  , xpd     = TRUE
                  , bty     = "n"
                  )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Open and set plot window.                                               #
            #------------------------------------------------------------------------------#
            par(mar=c(3.1,4.6,3.1,1.6))
            plot.new()
            plot.window(xlim=sm.xlimit,ylim=sm.ylimit)
            axis(side=1,las=1,at=em.pretty$levels,labels=em.pretty$labels)
            axis(side=2,las=1)
            title(main=le.cmean,ylab=ley,cex.main=1.0)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #       Plot range bands.                                                      #
            #------------------------------------------------------------------------------#
            polygon( x       = c( tomonth[eshow]
                                , rev(tomonth[eshow])
                                , NA
                                , tomonth[eshow]
                                , rev(tomonth[eshow])
                                )#end c
                   , y       = c( cmean[eshow,"qlwr","test"]
                                , rev(cmean[eshow,"qupr","test"])
                                , NA
                                , cmean[eshow,"qlwr","ctrl"]
                                , rev(cmean[eshow,"qupr","ctrl"])
                                )#end c
                   , col     = cband.col
                   , angle   = cband.angle
                   , density = cband.density
                   , lwd     = cband.lwd
                   )#end polygon
            #------------------------------------------------------------------------------#



            #------ Plot monthly averages. ------------------------------------------------#
            lines( x    = tomonth[eshow]
                 , y    = cmean[eshow,"mean","test"]
                 , type = "l"
                 , col  = cmean.col["test"]
                 , lty  = cmean.lty["test"]
                 , lwd  = 2.5
                 )#end lines
            lines( x    = tomonth[eshow]
                 , y    = cmean[eshow,"mean","ctrl"]
                 , type = "l"
                 , col  = cmean.col["ctrl"]
                 , lty  = cmean.lty["ctrl"]
                 , lwd  = 2.5
                 )#end lines
            box()
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #     Plot vertical profile.                                                         #
      #------------------------------------------------------------------------------------#
         mon.loop = sequence(12*this.hgtprof)
         for (m in mon.loop){
            #----- Aliases for current month. ---------------------------------------------#
            mlabel   = paste0("M",sprintf("%2.2i",m))
            mdesc    = month.name[m]
            mzsuffix = paste0(this.vnam,"-",iata,"_",mlabel)
            le.mprof = paste0(longname,"\n","Vertical profile (",mdesc,") ")
            le.xprof = desc.unit(desc=this.desc,unit=this.unit)
            le.yprof = desc.unit(desc="Height" ,unit=untab$m  )
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot monthly mean diurnal cycle.                                        #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = file.path(outhpmain,paste0("hgtprof-",mzsuffix,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = xsize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               layout(mat= rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #----- Plot patch properties. ----------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = cband.leg
                     , fill    = cband.col
                     , border  = cband.col
                     , angle   = cband.angle
                     , density = cband.density
                     , col     = cmean.col
                     , lty     = cmean.lty
                     , lwd     = 2.5
                     , ncol    = 2
                     , cex     = 1.0
                     , xpd     = TRUE
                     , bty     = "n"
                     )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Open device.                                                         #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.6,3.1,1.6))
               plot.new()
               plot.window(xlim=hm.xlimit,ylim=hm.ylimit)
               axis(side=1,las=1,at=hm.xpretty,labels=hm.xlabels)
               axis(side=2,las=1,at=hm.ypretty,labels=hm.ylabels)
               title(main=le.mprof,xlab=le.xprof,ylab=le.yprof,cex.main=1.0)
               #---------------------------------------------------------------------------#



               #------ Plot range bands. --------------------------------------------------#
               polygon( x       = c(hprof[m,,"qlwr","test"],rev(hprof[m,,"qupr","test"]))
                      , y       = c(hgt.lyr.levs,rev(hgt.lyr.levs))
                      , col     = cband.col    ["test"]
                      , angle   = cband.angle  ["test"]
                      , density = cband.density["test"]
                      , lwd     = cband.lwd    ["test"]
                      )#end polygon
               polygon( x       = c(hprof[m,,"qlwr","ctrl"],rev(hprof[m,,"qupr","ctrl"]))
                      , y       = c(hgt.lyr.levs,rev(hgt.lyr.levs))
                      , col     = cband.col    ["ctrl"]
                      , angle   = cband.angle  ["ctrl"]
                      , density = cband.density["ctrl"]
                      , lwd     = cband.lwd    ["ctrl"]
                      )#end polygon
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Plot vertical profile.                                               #
               #---------------------------------------------------------------------------#
               lines( x    = hprof[m,,"mean","test"]
                    , y    = hgt.lyr.levs
                    , type = "l"
                    , col  = cmean.col["test"]
                    , lty  = cmean.lty["test"]
                    , lwd  = 2.5
                    )#end lines
               lines( x    = hprof[m,,"mean","ctrl"]
                    , y    = hgt.lyr.levs
                    , type = "l"
                    , col  = cmean.col["ctrl"]
                    , lty  = cmean.lty["ctrl"]
                    , lwd  = 2.5
                    )#end lines
               box()
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end for (m in mon.loop)
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #     Plot vertical profile.                                                         #
      #------------------------------------------------------------------------------------#
         mon.loop = sequence(12*this.hgtprof)
         for (m in mon.loop){
            #----- Aliases for current month. ---------------------------------------------#
            mlabel   = paste0("M",sprintf("%2.2i",m))
            mdesc    = month.name[m]
            mzsuffix = paste0(this.vnam,"-",iata,"_",mlabel)
            le.mprof = paste0(longname,"\n","Vertical profile (",mdesc,") ")
            le.xprof = desc.unit(desc=this.desc,unit=this.unit)
            le.yprof = desc.unit(desc="Height" ,unit=untab$m  )
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot monthly mean diurnal cycle.                                        #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = file.path(outhcmain,paste0("cumprof-",mzsuffix,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = xsize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               layout(mat= rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #----- Plot patch properties. ----------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = cband.leg
                     , fill    = cband.col
                     , border  = cband.col
                     , angle   = cband.angle
                     , density = cband.density
                     , col     = cmean.col
                     , lty     = cmean.lty
                     , lwd     = 2.5
                     , ncol    = 2
                     , cex     = 1.0
                     , xpd     = TRUE
                     , bty     = "n"
                     )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Open device.                                                         #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.6,3.1,1.6))
               plot.new()
               plot.window(xlim=cm.xlimit,ylim=cm.ylimit)
               axis(side=1,las=1,at=cm.xpretty,labels=cm.xlabels)
               axis(side=2,las=1,at=cm.ypretty,labels=cm.ylabels)
               title(main=le.mprof,xlab=le.xprof,ylab=le.yprof,cex.main=1.0)
               #---------------------------------------------------------------------------#



               #------ Plot range bands. --------------------------------------------------#
               polygon( x       = c(cprof[m,,"qlwr","test"],rev(cprof[m,,"qupr","test"]))
                      , y       = c(hgt.cdf.levs,rev(hgt.cdf.levs))
                      , col     = cband.col    ["test"]
                      , angle   = cband.angle  ["test"]
                      , density = cband.density["test"]
                      , lwd     = cband.lwd    ["test"]
                      )#end polygon
               polygon( x       = c(cprof[m,,"qlwr","ctrl"],rev(cprof[m,,"qupr","ctrl"]))
                      , y       = c(hgt.cdf.levs,rev(hgt.cdf.levs))
                      , col     = cband.col    ["ctrl"]
                      , angle   = cband.angle  ["ctrl"]
                      , density = cband.density["ctrl"]
                      , lwd     = cband.lwd    ["ctrl"]
                      )#end polygon
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Plot vertical profile.                                               #
               #---------------------------------------------------------------------------#
               lines( x    = cprof[m,,"mean","test"]
                    , y    = hgt.cdf.levs
                    , type = "l"
                    , col  = cmean.col["test"]
                    , lty  = cmean.lty["test"]
                    , lwd  = 2.5
                    )#end lines
               lines( x    = cprof[m,,"mean","ctrl"]
                    , y    = hgt.cdf.levs
                    , type = "l"
                    , col  = cmean.col["ctrl"]
                    , lty  = cmean.lty["ctrl"]
                    , lwd  = 2.5
                    )#end lines
               box()
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end for (m in mon.loop)
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#









      #------------------------------------------------------------------------------------#
      #     Plot maps using emean.                                                         #
      #------------------------------------------------------------------------------------#
      echron = chron(dimnames(emean)[[1]])
      ey4m2  = paste0(sprintf("%4.4i",numyears(echron))
                     ,"-"
                     ,sprintf("%2.2i",nummonths(echron))
                     )#end paste0
      emzy4  = paste(month.name[nummonths(echron)],sprintf("%4.4i",numyears(echron)))
      for (e in which(eshow)){
         #----- Alias for current time suffix (to be appended to file names). -------------#
         esuffix = paste0(this.vnam,"-",iata,"_",ey4m2[e])
         etitle  = paste0(this.desc," - ",emzy4[e])
         #---------------------------------------------------------------------------------#


         #------ Get indices for maps. ----------------------------------------------------#
         ic = cbind(e,ptable$ipa[model$imap],5)
         it = cbind(e,ptable$ipa[model$imap],ptable$ilight[model$imap])
         #---------------------------------------------------------------------------------#


         #------ Make maps for control and test. ------------------------------------------#
         ectrl = matrix(data=emean[ic],nrow=n.xy,ncol=n.xy)
         etest = matrix(data=emean[it],nrow=n.xy,ncol=n.xy)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Check that this month has at least one valid (non-NA) point.                #
         #---------------------------------------------------------------------------------#
         if (all(is.finite(mm.yrange[,e]))){
            #------------------------------------------------------------------------------#
            #     Prepare settings for maps.                                               #
            #------------------------------------------------------------------------------#
            mm.zat    = pretty(x=mm.yrange[,e])
            mm.zlab   = sprintf("%g",mm.zat)
            mm.zlimit = range(mm.zat)
            mm.zlevs  = seq(from=mm.zlimit[1],to=mm.zlimit[2],length.out=ncpbks)
            mm.zcols  = this.cscheme(n=ncolpal)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot monthly means.                                                     #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = file.path(outmpmain,paste0("map-",esuffix,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = msize
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout(mat= rbind(c(2,3,1)),widths=c(rep(x=(2.-f.leg)/4,times=2),f.leg/2))
               #---------------------------------------------------------------------------#



               #----- Plot colour palette. ------------------------------------------------#
               par(mar=c(1.1,0.5,2.6,3.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=mm.zlimit,xaxs="i",yaxs="i")
               rect( xleft   = 0
                   , ybottom = mm.zlevs[-ncpbks]
                   , xright  = 1
                   , ytop    = mm.zlevs[-1]
                   , col     = mm.zcols
                   , border  = mm.zcols
                   )#end rect
               axis(side=4,las=1,at=mm.zat,labels=mm.zlab)
               title( main     = desc.unit(desc=NULL,unit=untab[[this.unit]])
                    , cex.main = 1.0
                    , line     = 1.0
                    )#end title
               #---------------------------------------------------------------------------#



               #----- Plot homogeneous illumination. --------------------------------------#
               par(mar=c(1.1,1.1,2.6,0.1))
               plot.new()
               plot.window(xlim=xy.range,ylim=xy.range,xaxs="i",yaxs="i")
               image( x         = xy.gaps
                    , y         = xy.gaps
                    , z         = ectrl
                    , col       = mm.zcols
                    , add       = TRUE
                    , breaks    = mm.zlevs
                    , useRaster = TRUE
                    )#end image
               box()
               title(main="Homogeneous illumination",cex.main=1.25,line=0.5)
               #---------------------------------------------------------------------------#



               #----- Plot heterogeneous illumination. ------------------------------------#
               par(mar=c(1.1,1.1,2.6,0.1))
               plot.new()
               plot.window(xlim=xy.range,ylim=xy.range,xaxs="i",yaxs="i")
               image( x         = xy.gaps
                    , y         = xy.gaps
                    , z         = etest
                    , col       = mm.zcols
                    , add       = TRUE
                    , breaks    = mm.zlevs
                    , useRaster = TRUE
                    )#end image
               box()
               title(main="Heterogeneous illumination",cex.main=1.25,line=0.5)
               #---------------------------------------------------------------------------#


               #----- Main title. ---------------------------------------------------------#
               mtext(text=etitle,side=3,outer=TRUE,cex=1.25,font=2)
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end if (all(is.finite(mm.yrange[,e])))
         #---------------------------------------------------------------------------------#
      }#end for (e in which(eshow))
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end for (s in sequence(nsites))
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
