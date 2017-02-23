#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
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
outroot = file.path(here,paste0("patch_comp_ibg",sprintf("%2.2i",ibackground)))
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.hour  = c(FALSE,TRUE)[1]
reload.range = c(FALSE,TRUE)[1]
rdata.path   = file.path(here,"RData_patches")
rdata.suffix = "patches_ed22.RData"
#------------------------------------------------------------------------------------------#



#----- Default settings. ------------------------------------------------------------------#
emean.yeara = 2008  # First year
emean.yearz = 2010  # Last year
eshow.yeara = 2008  # First year to show
eshow.yearz = 2008  # First year to show
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Quantiles that determine which patches to show.                                       #
#------------------------------------------------------------------------------------------#
qshow = c(0,0.25,0.50,0.75,1.00)
#------------------------------------------------------------------------------------------#



#----- Map settings. ----------------------------------------------------------------------#
pat.min.area = 0.0025
n.min.area   = 20
dxy.gap      = 20
#------------------------------------------------------------------------------------------#



#----- Confidence band. ----------------------------------------------------------------------#
cband         = -1 # 2*pnorm(q=1)-1
cband.leg     = c(test="Heterogeneous",ctrl="Homogeneous")
cband.col     = c(test="#1E90FF"      ,ctrl="#BFBFBF"    )
cmean.col     = c(test="#191970"      ,ctrl="#404040"    )
cband.angle   = c(test=-45            ,ctrl=+45          )
cband.density = c(test= 25            ,ctrl= 25          )
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
#              Logical   -- If TRUE, it uses all sites; FALSE doesn't make sense...        #
#              Character -- vector with the 3-letter code ("iata") of selected sites       #
#              Numeric   -- Site index.                                                    #
#                                                                                          #
#------------------------------------------------------------------------------------------#
eort       = "t"
n          = 0
sites      = list()
n          = n + 1
sites[[n]] = list( iata = "tnf_iage060"
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
compvar[[ n]] = list( vnam     = "gpp"
                    , desc     = "Gross primary productivity"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "npp"
                    , desc     = "Net primary productivity"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "plant.resp"
                    , desc     = "Plant respiration"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "iclife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "cba"
                    , desc     = "Carbon balance"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "clife"
                    , qmean    = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "reco"
                    , desc     = "Ecosystem respiration"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "iclife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "nep"
                    , desc     = "Net Ecosystem Productivity"
                    , unit     = "kgcom2oyr"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "hflxca"
                    , desc     = "Sensible heat flux"
                    , unit     = "wom2"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "wflxca"
                    , desc     = "Water vapour flux"
                    , unit     = "kgwom2oday"
                    , cscheme  = "ipanoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "transp"
                    , desc     = "Transpiration"
                    , unit     = "kgwom2oday"
                    , cscheme  = "ipanoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "rshortup"
                    , desc     = "Upward SW radiation"
                    , unit     = "wom2"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "rlongup"
                    , desc     = "Upward LW radiation"
                    , unit     = "wom2"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "parup"
                    , desc     = "Upward PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.gnd"
                    , desc     = "Ground absorption - PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "rshort.gnd"
                    , desc     = "Ground absorption - SW"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "sm.stress"
                    , desc     = "Soil moisture stress"
                    , unit     = "empty"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.gpp"
                    , desc     = "Leaf GPP"
                    , unit     = "kgcom2loyr"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.temp"
                    , desc     = "Mean Leaf Temperature"
                    , unit     = "degC"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.vpd"
                    , desc     = "Mean Leaf VPD"
                    , unit     = "hpa"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.gsw"
                    , desc     = "Stomatal conductance"
                    , unit     = "kgwom2loday"
                    , cscheme  = "ipanoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.leaf"
                    , desc     = "Leaf Absorption - PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.leaf.beam"
                    , desc     = "Leaf Absorption - Direct PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "par.leaf.diff"
                    , desc     = "Leaf Absorption - Diffuse PAR"
                    , unit     = "umolom2os"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.par"
                    , desc     = "Norm. Leaf Absorption - PAR"
                    , unit     = "umolom2los"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.par.beam"
                    , desc     = "Norm. Leaf Absorption - Direct PAR"
                    , unit     = "umolom2los"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "leaf.par.diff"
                    , desc     = "Norm. Leaf Absorption - Diffuse PAR"
                    , unit     = "umolom2los"
                    , cscheme  = "panoply"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "assim.light"
                    , desc     = "Light-limited Assimilation"
                    , unit     = "umolom2los"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "assim.rubp"
                    , desc     = "RuBP-limited Assimilation"
                    , unit     = "umolom2los"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "assim.co2"
                    , desc     = "CO2-limited Assimilation"
                    , unit     = "umolom2los"
                    , cscheme  = "clife"
                    , qmean    = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.lpar"
                    , desc     = "Daytime PAR absorportion by leaves"
                    , unit     = "hpa"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.lvpd"
                    , desc     = "Daytime Leaf VPD"
                    , unit     = "hpa"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.ltemp"
                    , desc     = "Daytime Leaf Temperature"
                    , unit     = "degC"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.sms"
                    , desc     = "Daytime soil moisture stress"
                    , unit     = "empty"
                    , cscheme  = "panoply"
                    , qmean    = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam     = "phap.lgsw"
                    , desc     = "Daytime stomatal conductance"
                    , unit     = "kgwom2loday"
                    , cscheme  = "ipanoply"
                    , qmean    = FALSE
                    )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.suffix  = "imetrad01"
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

depth             = 300                  # PNG resolution, in pixels per inch
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
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
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



#------ Find out whether to plot annual means (it must have at least 2 years... -----------#
plot.ymean = eshow.yearz > eshow.yeara
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
outsimul = file.path(outroot,sim.suffix)
outymean = file.path(outsimul,"tseries_ymean")
outemean = file.path(outsimul,"tseries_emean")
outcmean = file.path(outsimul,"tseries_cmean")
outqmean = file.path(outsimul,"tseries_qmean")
outmmean = file.path(outsimul,"maps_emean")
if (! file.exists(outroot )             ) dir.create(outroot )
if (! file.exists(outsimul)             ) dir.create(outsimul)
if (! file.exists(outymean) & plot.ymean) dir.create(outymean)
if (! file.exists(outemean)             ) dir.create(outemean)
if (! file.exists(outcmean)             ) dir.create(outcmean)
if (! file.exists(outqmean)             ) dir.create(outqmean)
if (! file.exists(outmmean)             ) dir.create(outmmean)
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
      cat0(" + Load data from ",basename(rdata.iata),"...")
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
if (length(loop.sites) != 0) cat0(" + Processing missing hourly data...")
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
   cat0("   - Site :",this$longname,"...")
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Get all the statistics and actual values for every simulation.                    #
   #---------------------------------------------------------------------------------------#
   cat0("    * Aggregate and find statistics for simulations for this site...")

   #----- Load hourly averages (Control). -------------------------------------------------#
   ctrl.name = paste0(eort,iata,"_",sim.struct["ctrl"],"_",sim.suffix)
   cat0("      # Control: ",ctrl.name,"...")
   ctrl.path = file.path(here,ctrl.name)
   ctrl.file = file.path(ctrl.path,"rdata_month",paste0(ctrl.name,".RData"))
   load(ctrl.file)
   ctrl      = datum
   rm(datum)
   #---------------------------------------------------------------------------------------#


   #----- Load hourly averages. -----------------------------------------------------------#
   test.name = paste0(eort,iata,"_",sim.struct["test"],"_",sim.suffix)
   cat0("      # Test: ",test.name,"...")
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
      cat0("      # Test table: ",basename(test.name),"...")
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
   qfbeam          = quantile(x=wfbeam,probs=qshow,na.rm=TRUE)
   qidx.show       = mapply(FUN=which.closest,x=qfbeam,MoreArgs=list(A=wfbeam))
   model$qidx.show = qidx.show
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Define matrix placeholder.                                                        #
   #---------------------------------------------------------------------------------------#
   e.empty = array( data     = NA
                  , dim      = c(nemean,npatches,nlight)
                  , dimnames = list( as.character(model$tomonth)
                                   , sequence(npatches)
                                   , light.key
                                   )#end list
                  )#end array
   c.empty = array( data     = NA
                  , dim      = c(nemean,3,2)
                  , dimnames = list( as.character(model$tomonth)
                                   , c("qlwr","mean","qupr")
                                   , c("ctrl","test")
                                   )#end list
                  )#end array
   q.empty = array( data     = NA
                  , dim      = c(nemean,ndcycle,npatches,nlight)
                  , dimnames = list( as.character(model$tomonth)
                                   , dclabel
                                   , sequence(npatches)
                                   , light.key
                                   )#end list
                  )#end array
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Load all variables, and append them to matrices.                                 #
   #---------------------------------------------------------------------------------------#
   cat0("       ~ Load variables...")
   for (v in sequence(ncompvar)){
      #----- Load information. ------------------------------------------------------------#
      this.vnam      = compvar$vnam [v]
      this.desc      = compvar$desc [v]
      this.unit      = compvar$unit [v]
      this.qmean     = compvar$qmean[v]
      cat0("         > ",this.desc,"...")
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Initialise array for this variable.                                            #
      #------------------------------------------------------------------------------------#
      model[[this.vnam]] = list( emean = e.empty, qmean = q.empty, cmean = c.empty)
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
               csdev = weighted.sd      (x=cpatch[[stamp]],w=carea[[stamp]]        )
               model[[this.vnam]]$cmean[e,,"ctrl"] = cmean + c(-1.,0.,1.)*abs(cband)*csdev
            }#end if (cband > 0)
            #------------------------------------------------------------------------------#
         }#end if (! is.null(cpatch))
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
               tsdev = weighted.sd      (x=tpatch[[stamp]],w=tarea[[stamp]]        )
               model[[this.vnam]]$cmean[e,,"test"] = tmean + c(-1.,0.,1.)*abs(cband)*tsdev
            }#end if (cband > 0)
            #------------------------------------------------------------------------------#
         }#end if (! is.null(tpatch)){
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
         }#end if (! is.null(tpatch)){
         #---------------------------------------------------------------------------------#
      }#end for (e in emean.loop)
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
   cat0(" + Saving processed data to ",basename(rdata.iata),"...")
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
cat0(" + Plot long-term time series...")
for (s in sequence(nsites)){
   #----- Loop over variables. ------------------------------------------------------------#
   iata     = sites$iata[s]
   longname = sites$desc[s]
   model    = res[[iata]]
   npatches = model$npatches
   props    = model$props
   ptable   = model$ptable
   #---------------------------------------------------------------------------------------#
   
   cat0("   - ",longname,"...")

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
      tomonth         = model$tomonth
      toyear          = sort(unique(model$toyear))
      emean           = model[[this.vnam]]$emean
      cmean           = model[[this.vnam]]$cmean
      ymean           = qapply(X=emean,INDEX=model$toyear,DIM=1,FUN=mean,na.rm=TRUE)
      qmean           = model[[this.vnam]]$qmean
      wshow           = model$tomonth %wr% when.show
      qidx.show       = model$qidx.show
      cat0("     > ",this.desc,"...")
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


      #------------------------------------------------------------------------------------#
      #     Find plot limits for time.                                                     #
      #------------------------------------------------------------------------------------#      
      em.xlimit = chron(range(x=model$tomonth[wshow]))
      cm.xlimit = em.xlimit
      ym.xlimit = pretty.xylim(u=model$toyear [wshow])
      qm.xlimit = pretty.xylim(u=model$dchour)
      em.pretty = pretty.time (when=em.xlimit)
      cm.pretty = em.pretty
      qm.pretty = pretty.elapsed(x=model$dchour,base=24)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find range for each patch.                                                     #
      #------------------------------------------------------------------------------------#
      eshow     = chron(dimnames(emean)[[1]])      %wr% when.show
      cshow     = chron(dimnames(cmean)[[1]])      %wr% when.show
      yshow     = as.numeric(dimnames(ymean)[[1]]) %wr% numyears(when.show)
      em.yrange = apply(X=emean[eshow,,,drop=FALSE],MARGIN=2,FUN=range,finite=TRUE)
      cm.yrange = range(c(cmean[eshow,,]),finite=TRUE)
      ym.yrange = apply(X=ymean[yshow,,,drop=FALSE],MARGIN=2,FUN=range,finite=TRUE)
      qm.yrange = apply(X=qmean                    ,MARGIN=3,FUN=range,finite=TRUE)
      mm.yrange = apply(X=emean                    ,MARGIN=1,FUN=range,finite=TRUE)
      #------------------------------------------------------------------------------------#




      #------ Set some common features. ---------------------------------------------------#
      ley     = desc.unit(desc=this.desc,unit=this.unit)
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
         outmpmain = file.path(outmmean,this.vnam)
         if (! file.exists(outmpmain)) dir.create(outmpmain)
         #---------------------------------------------------------------------------------#


         #----- Limits for this patch. ----------------------------------------------------#
         em.ylimit = pretty.xylim(u=c(em.yrange[,p]))
         ym.ylimit = pretty.xylim(u=c(ym.yrange[,p]))
         qm.ylimit = pretty.xylim(u=c(qm.yrange[,p]))
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
            if (outform[o] %in% "x11"){
               X11(width=xsize$width,height=xsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=xsize$width,height=xsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=xsize$width,height=xsize$height
                         ,pointsize=ptsz,paper=xsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=xsize$width,height=xsize$height
                  ,pointsize=ptsz,paper=xsize$paper)
            }#end if
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
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
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
            if (outform[o] %in% "x11"){
               X11(width=xsize$width,height=xsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=xsize$width,height=xsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=xsize$width,height=xsize$height
                         ,pointsize=ptsz,paper=xsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=xsize$width,height=xsize$height
                  ,pointsize=ptsz,paper=xsize$paper)
            }#end if
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
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
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
               if (outform[o] %in% "x11"){
                  X11(width=xsize$width,height=xsize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=xsize$width,height=xsize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=xsize$width,height=xsize$height
                            ,pointsize=ptsz,paper=xsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=xsize$width,height=xsize$height
                     ,pointsize=ptsz,paper=xsize$paper)
               }#end if
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
               if (outform[o] %in% c("x11","quartz")){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               dummy = clean.tmp()
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
         cm.ylimit = pretty.xylim(u=c(cm.yrange    ))
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
            if (outform[o] %in% "x11"){
               X11(width=xsize$width,height=xsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=xsize$width,height=xsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=xsize$width*depth,height=xsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=xsize$width,height=xsize$height
                         ,pointsize=ptsz,paper=xsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=xsize$width,height=xsize$height
                  ,pointsize=ptsz,paper=xsize$paper)
            }#end if
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
                  , lty     = "solid"
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
            plot.window(xlim=cm.xlimit,ylim=cm.ylimit)
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
                   )#end polygon
            #------------------------------------------------------------------------------#



            #------ Plot monthly averages. ------------------------------------------------#
            lines( x    = tomonth[eshow]
                 , y    = cmean[eshow,"mean","test"]
                 , type = "l"
                 , col  = cmean.col["test"]
                 , lwd  = 2.5
                 )#end lines
            lines( x    = tomonth[eshow]
                 , y    = cmean[eshow,"mean","ctrl"]
                 , type = "l"
                 , col  = cmean.col["ctrl"]
                 , lwd  = 2.5
                 )#end lines
            box()
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
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
         #     Prepare settings for maps.                                                  #
         #---------------------------------------------------------------------------------#
         mm.zat    = pretty(x=mm.yrange[,e])
         mm.zlab   = sprintf("%g",mm.zat)
         mm.zlimit = range(mm.zat)
         mm.zlevs  = seq(from=mm.zlimit[1],to=mm.zlimit[2],length.out=ncpbks)
         mm.zcols  = this.cscheme(n=ncolpal)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Plot monthly means.                                                        #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            fichier = file.path(outmpmain,paste0("map-",esuffix,".",outform[o]))
            if (outform[o] %in% "x11"){
               X11(width=msize$width,height=msize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=msize$width,height=msize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=msize$width*depth,height=msize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=msize$width*depth,height=msize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=msize$width,height=msize$height
                         ,pointsize=ptsz,paper=msize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=msize$width,height=msize$height
                  ,pointsize=ptsz,paper=msize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2.5,0))
            layout(mat= rbind(c(2,3,1)),widths=c(rep(x=(2.-f.leg)/4,times=2),f.leg/2))
            #------------------------------------------------------------------------------#



            #----- Plot colour palette. ---------------------------------------------------#
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
            title(main=desc.unit(desc=NULL,unit=untab[[this.unit]]),cex.main=1.0,line=1.0)
            #------------------------------------------------------------------------------#



            #----- Plot homogeneous illumination. -----------------------------------------#
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
            #------------------------------------------------------------------------------#



            #----- Plot heterogeneous illumination. ---------------------------------------#
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
            #------------------------------------------------------------------------------#


            #----- Main title. ------------------------------------------------------------#
            mtext(text=etitle,side=3,outer=TRUE,cex=1.25,font=2)
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] %in% c("x11","quartz")){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
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
