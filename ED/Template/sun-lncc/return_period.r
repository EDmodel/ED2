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
here           = getwd()                          # Current directory
srcdir         = "/prj/prjidfca/marcosl/Util/Rsc" # Script directory
stext.default  = "stext16"                   # Default soil texture
drain.default  = "r+000"                     # Default rainfall
ibackground    = 0                           # Target background colour:
                                             #   (to adjust foreground colours accordingly)
                                             # 0 -- White
                                             # 1 -- Pitch black
                                             # 2 -- Dark grey
bg.default     = paste("ibg",sprintf("%2.2i",ibackground),sep="")
outroot        = file.path(here,paste("vulnerability",bg.default,sep="_"))
comp.prefix    = "stext"
crit.loss      = 0.2                         # Critical loss of biomass (% of initial)
n.sample       = 1000                        # Number of samples
n.pred         = 20000                       # Number of points for curve drawing
std.pret.limit = c(1,500)                    # Fixed return period
nls.optim      = TRUE                        # Use nls for optimisation (FALSE uses optim)
robust         = TRUE                        # Use robust fitting in the NLS
skew.optim     = FALSE                       # Use skew normal distribution for residuals?
                                             #     (FALSE uses normal distribution aka 
                                             #      least squares maximisation)
verbose.optim  = FALSE                       # Dump optimisation progress on screen?
                                             # 0 or FALSE -- No
                                             # 1 or TRUE  -- Optimiser main step
                                             # 2 or > 2   -- Optimiser step and sub-step
save.every     = 2880                        # Save the partially loaded data every...
pch.return.type = FALSE                      # PCH to denote the return period type?
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
rdata.path       = file.path(here,"RData_scenario") # Path with R object.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Variable name.template has the general format of the simulations.  Put some patterns #
# instead of the actual values, and use these patterns in the variables below to replace   #
# and substitute by the actual names.  For the other variables you must provide the        #
# following lists:                                                                         #
#                                                                                          #
# global   -- These are global variables, forced to be unidimensional.  Direct comparison  #
#             plots will never happen amongst these variables.  At least one global        #
#             dimension must be given.                                                     #
# panel    -- Similar to global variables, they will be used to generate multiple panels   #
#             0  dimension  -- Single panel                                                #
#             1  dimension  -- One panel for each value                                    #
#             2  dimensions -- All combinations                                            #
#             3+ dimension  -- Not allowed                                                 #
# scenario -- Scenario variables.                                                          #
#             0  dimension  -- Not allowed (and btw, did you notice the script is called   #
#                              compare scenarios?)                                         #
#             1  dimension  -- line plots with the average properties                      #
#             2  dimensions -- image plots.                                                #
#             3+ dimensions -- Not allowed. I'm not the one who likes 3D plots :).         #
#                                                                                          #
#    Each list element of the list must contain the following variables:                   #
#                                                                                          #
#    key     -- the identifier for the variable, that is in the name of the simulation     #
#               and may be used in the output file name and in the array dimension         #
#               names.                                                                     #
#    desc    -- a nicer name that describes the variable.  They will be used in title      #
#               or axis labels.                                                            #
#    pattern -- the pattern that will be replaced by the key in name.template.  This       #
#               should be a scalar.                                                        #
#                                                                                          #
#    The following elements are used only by the scenario variables, but they must be      #
#       set for all elements (NA is enough for the others)                                 #
#                                                                                          #
#    value   -- The values so we can build the image plot with reasonable axis.            #
#    default -- The default value, the other simulations will be scaled against this       #
#               variable.                                                                  #
#    colour  -- The colours that identify this scenario.                                   #
#    pch     -- The symbols that identify this scenario.                                   #
#    alabel  -- The axis labels.  This should be a scalar.                                 #
#                                                                                          #
#                                                                                          #
#                                                                                          #
#    In addition, you must define a list with the configurations of each realisation.      #
# In case you didn't run any realisation, define the elements as empty vectors.            #
# (only key, desc, and pattern are defined).                                               #
#------------------------------------------------------------------------------------------#
name.template    = "tPPP_rRRRR_tTTTT_real-ZZ_iphenDDD_stextSS"
use.global       = 1   # Which global to use (TRUE means all of them)
#----- Global variables. ------------------------------------------------------------------#
global         = list()
#global$stext   = list( key     = c("stext02","stext06","stext08","stext16","stext11")
#                     , desc    = paste( "Soil type:"
#                                      , c("Loamy sand","Sandy clay loam","Clayey loam"
#                                         ,"Clayey sand","Clay")
#                                      ,sep = " "
#                                      )#end paste
#                     , pattern = "stextSS"
#                     , value   = NA_real_
#                     , label   = NA_real_
#                     , default = NA_character_
#                     , colour  = NA_character_
#                     , pch     = NA_integer_
#                     , alabel  = NA_character_
#                     )#end list
global$dtemp   = list( key     = c("t+000","t+100","t+200","t+300")
                     , desc    = c("","","","")
                     , legend  = c("","","","")
                     , parse   = FALSE
                     , pattern = "tTTTT"
                     , default = NA_integer_
                     , value   = NA_real_
                     , label   = NA_real_
                     , colour  = NA_character_
                     , pch     = NA_integer_
                     , alabel  = NA_character_
                     )#end list
#----- Panel variables. -------------------------------------------------------------------#
panel          = list()
panel$iata     = list( key     = c("gyf","s67")
                     , desc    = c("Paracou","Santarem")
                     , legend  = c("Paracou","Santarem")
                     , parse   = FALSE
                     , pattern = "PPP"
                     , default = NA_character_
                     , value   = NA_real_
                     , label   = NA_real_
                     , colour  = NA_character_
                     , pch     = NA_integer_
                     , alabel  = NA_character_
                     )#end list
panel$iphen    = list( key     = c("iphen-01","iphen+02")
                     , desc    = paste("Phenology:",c("Evergreen","Drought Deciduous"))
                     , legend  = paste("Phenology:",c("Evergreen","Drought Deciduous"))
                     , parse   = FALSE
                     , pattern = "iphenDDD"
                     , default = NA_character_
                     , value   = NA_real_
                     , label   = NA_real_
                     , colour  = NA_character_
                     , pch     = NA_integer_
                     , alabel  = NA_character_
                     )#end list
#----- Scenario variables. ----------------------------------------------------------------#
scenario       = list()
scenario$drain = list( key     = c("r+000","r-020","r-040","r-060","r-080","r-100"
                                  ,"r-120","r-140","r-160")
                     , desc    = c("dR =  0.0S","dR = -0.2S","dR = -0.4S","dR = -0.6S"
                                  ,"dR = -0.8S","dR = -1.0S","dR = -1.2S","dR = -1.4S"
                                  ,"dR = -1.6S")
                     , legend  = c("Delta*xi ==  0.0*omega","Delta*xi == -0.2*omega"
                                  ,"Delta*xi == -0.4*omega","Delta*xi == -0.6*omega"
                                  ,"Delta*xi == -0.8*omega","Delta*xi == -1.0*omega"
                                  ,"Delta*xi == -1.2*omega","Delta*xi == -1.4*omega"
                                  ,"Delta*xi == -1.6*omega")
                     , parse   = TRUE
                     , pattern = "rRRRR"
                     , default = drain.default
                     , value   = c(  0.0,  0.2,  0.4,  0.6,  0.8,  1.0,  1.2,  1.4,  1.6)
                     , label   = c(  0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6)
                     , colour  = c("#003264","#0082C8","#46B4FF","#B4E6FF","#DAF7F1"
                                  ,"#E6E6B4","#FFB43C","#C85A0A","#960000")
                     , pch     = c(12L,9L,10L,0L,5L,1L,2L,6L,8L)
                     , alabel  = c("Mean rainfall change [Scale]")
                    )#end list
scenario$stext = list( key     = c("stext02","stext06","stext08","stext16","stext11")
                     , desc    = c("Loamy sand","Sandy clay loam","Clay loam"
                                  ,"Clayey sand","Clay")
                     , legend  = c("Loamy sand","Sandy clay loam","Clay loam"
                                  ,"Clayey sand","Clay")
                     , parse   = FALSE
                     , pattern = "stextSS"
                     , default = stext.default
                     , value   = c(1,2,3,4,5)
                     , label   = c("LSa","SaCL","CL","CSa","C")
                     , colour  = c("#3B24B3","#2996CC","#A3CC52","#E65C17","#990F0F")
                     , pch     = c(12L,13L,5L,6L,8L)
                     , alabel  = c("Soil texture")
                     )#end list
#----- Realisation variables. -------------------------------------------------------------#
realisation = list( key     = paste("real",sprintf("%2.2i",0:15),sep="-")
                  , desc    = paste("Realisation",sprintf("%2.2i",0:15),sep=" ")
                  , pattern = "real-ZZ"
                  )#end list
#------------------------------------------------------------------------------------------#





#------ Miscellaneous settings. -----------------------------------------------------------#
yeara          = 1952         # First year we will include
yeare          = 1982         # First year to use in the averaged output
yearz          = 2011         # Last year we will include
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      = 3            # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                              # 3 -- 0-10; 10-35; 35-70; > 70 (cm)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform        = c("pdf")    # Formats for output file.  Supported formats are:
                             #   - "X11" - for printing on screen
                             #   - "eps" - for postscript printing
                             #   - "png" - for PNG printing
                             #   - "pdf" - for PDF printing

byeold         = TRUE        # Remove old files of the given format?

depth          = 96          # PNG resolution, in pixels per inch
paper          = "letter"    # Paper size, to define the plot shape
wpaper         = "long"      # Paper size for wide plots.
ptsz           = 16          # Font size.
lwidth         = 2.5         # Line width
plotgrid       = TRUE        # Should I plot the grid in the background? 

legwhere       = "topleft"   # Where should I place the legend?
inset          = 0.01        # Inset between legend and edge of plot region.
fracexp        = 0.40        # Expansion factor for y axis (to fit legend)
n.colourbar    = 32          # Number of colours for the colour bars
n.whitebar     =  1          # Number of levels around zero to be set to white
notch          = FALSE       # Add notches to the box plots.
mtext.xoff.im  = -9.00       # Offset for the x label
mtext.xoff     = -7.50       # Offset for the x label
mtext.xoff.e   = -4.50       # Offset for the x label
mtext.xoff.xyz = -3.90       # Offset for the x label
mtext.yoff     = -0.50       # Offset for the y label
mtext.yoff.im  = -1.50       # Offset for the y label
mtext.xadj     =  0.50       # Offset for the x label
mtext.yadj     =  0.65       # Offset for the y label
barplot.lwd    =  1.50       # Line width for the bar plots
plot.vec.axes  = FALSE       # Plot vector axes? (TRUE|FALSE)
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




#----- Load some packages and scripts. ----------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#
#      Function to find the default rainfall simulation associated with each run.          #
# Input variables:                                                                         #
#  x       -- the current simulation.  Only one simulation at a time.                      #
#  default -- the list of "default" rainfall simulations.                                  #
#------------------------------------------------------------------------------------------#
idx.zrain = function(x,default){
   #------ Sanity check on default. -------------------------------------------------------#
   if (is.data.frame(default)){
      z = default
   }else if(is.matrix(default)){
      z = as.data.frame(default)
   }else{
      stop("Default must be a matrix or a data frame!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #------ Sanity check on sizes. ---------------------------------------------------------#
   if (length(x) != ncol(z)){
      cat(" Size of x:               ",length(x),"\n")
      cat(" # of columns of default: ",ncol(z)  ,"\n")
      cat(" # of rows    of default: ",nrow(z)  ,"\n")
      stop(" Size of x must match number of columns of default!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #------ Make sure that g is a data.frame. ----------------------------------------------#
   g = as.data.frame( matrix( data     = x
                            , ncol     = ncol(z)
                            , nrow     = 1
                            , dimnames = list(NULL,dimnames(z)[[2]])
                            )#end matrix
                    )#end data.frame
   #---------------------------------------------------------------------------------------#



   #------ Find the correct index. --------------------------------------------------------#
   idx.iata  = which(z$iata        %in% g$iata )
   idx.iphen = which(z$iphen       %in% g$iphen)
   idx.stext = which(z$stext       %in% g$stext)
   idx       = intersect(idx.iata,intersect(idx.iphen,idx.stext))
   #---------------------------------------------------------------------------------------#
   #------ Crash in case idx doesn't have one element only. -------------------------------#
   if (length(idx) == 0){
      cat (" IDX.IATA  = ",idx.iata ,"\n")
      cat (" IDX.IPHEN = ",idx.iphen,"\n")
      cat (" IDX.STEXT = ",idx.stext,"\n")
      cat (" IDX       = ",idx      ,"\n")
      browser()
      stop(" IDX should have dimension 1!")
   }else{
      idx = idx[1]
   }#end if
   #---------------------------------------------------------------------------------------#

   return(idx)
}#end idx.zrain
#==========================================================================================#
#==========================================================================================#




#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#


#----- Make sure that the base directory exists. ------------------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size  = plotsize(proje=FALSE,paper=paper )
large = max(size$height,size$width)
ssize = modifyList(x=size,val=list(width=large,height=large,ratio=1))
psize = modifyList(x=size,val=list(height=size$height*3/2,ratio=size$ratio*2/3))
wsize = plotsize(proje=FALSE,paper=wpaper)
#------------------------------------------------------------------------------------------#



#----- Set some scaling for objects that may need to shrink for large font sizes. ---------#
cex.ptsz = min(15/ptsz,1.0)
#------------------------------------------------------------------------------------------#


#----- Set some dimensions associated with the simulations. -------------------------------#
n.global      = length(global         )
n.panel       = length(panel          )
n.scenario    = length(scenario       )
n.allscen     = n.scenario + 1
n.real        = length(realisation$key)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Check the number of dimensions in scenario and panel are correct.                    #
#------------------------------------------------------------------------------------------#
ok.global      = n.global > 0
ok.panel       = n.panel       %in% c(0,1,2)
ok.scenario    = n.scenario    %in% c(1,2)
if (! ( ok.global && ok.panel && ok.scenario)){
   cat (" Invalid dimensions!")
   cat (" Acceptable number of dimensions:")
   cat (" -- GLOBAL      = > 0"       ,"\n")
   cat (" -- PANEL       = 0, 1, or 2","\n")
   cat (" -- SCENARIO    = 1 or 2"    ,"\n")
   cat (" Given number of dimensions:")
   cat (" -- GLOBAL      =",n.global     ,"\n")
   cat (" -- PANEL       =",n.panel      ,"\n")
   cat (" -- SCENARIO    =",n.scenario   ,"\n")
   cat (" -- REALISATION =",n.real        ,"\n")
   stop("Correct the number of dimensions or re-write the script...")
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define some dimensions.                                                              #
#------------------------------------------------------------------------------------------#
pft.use       = c(   2,   3,   4,  18)          # PFTs to include (add PFT=18, the total)
pft.mp        = c(TRUE,TRUE,TRUE,TRUE)          # Include the PFT on multi-panel plots?
pft.key       = pft$key    [pft.use   ]         # PFT keys  (for dimnames)
pft.desc      = pft$name   [pft.use   ]         # PFT names (for titles)
pft.colour    = pft$colour [pft.use   ]         # PFT colours
dbh.use       = sequence(ndbh)                  # DBH classes to include
dbh.mp        = rep(TRUE,times=ndbh)            # DBH classes to plot on multi-panel
dbh.key       = dbhkeys    [dbh.use   ]         # DBH keys  (for dimnames)
dbh.desc      = dbhnames   [dbh.use   ]         # DBH names (for titles)
dbh.colour    = dbhcols    [dbh.use   ]         # DBH colours
season.use    = c(1,2,3,4,5)                    # Seasons to add (add season=5, the total)
season.key    = season.list[season.use]         # Keys for the seasons. 
season.desc   = season.full[season.use]         # Full names of all seasons.
season.colour = season.cols[season.use]         # Colours for seasons.
year.use      = yeara:yearz                     # Years to use
year.avg      = yeare:yearz                     # Years used for averages.
year.key      = year.use                        # Year keys  (for dimnames)
year.desc     = year.use                        # Year names (for titles)
#----- Same for drought length. -----------------------------------------------------------#
drought.use   = sequence(1)
drought.mp    = rep(TRUE,times=length(drought.use))
drought.key   = paste("dlen",sprintf("%2.2i",drought.use),sep="")
drought.desc  = paste(drought.use,"-yr drought",sep="")
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Transform pft.mp into an index linked to pft.use and associated variables.           #
#------------------------------------------------------------------------------------------#
pft.mp           = which(pft.mp)
pft.mp.key       = pft.key   [pft.mp]
pft.mp.desc      = pft.desc  [pft.mp]
pft.mp.colour    = pft.colour[pft.mp]
pft.bp           = sequence(length(pft.use)-1)
dbh.mp           = which(dbh.mp)
dbh.mp.key       = dbh.key   [dbh.mp]
dbh.mp.desc      = dbh.desc  [dbh.mp]
dbh.mp.colour    = dbh.colour[dbh.mp]
dbh.bp           = dbh.use
dbh.bp.key       = dbh.key   [dbh.bp]
dbh.bp.desc      = dbh.desc  [dbh.bp]
dbh.bp.colour    = dbh.colour[dbh.bp]
season.mp        = season.use   [-length(season.use)]
season.mp.key    = season.key   [-length(season.use)]
season.mp.desc   = season.desc  [-length(season.use)]
season.mp.colour = season.colour[-length(season.use)]
#------------------------------------------------------------------------------------------#


#------ Size of the useful dimensions. ----------------------------------------------------#
n.pft         = length(pft.use)          # Number of PFTs
n.pft.mp      = length(pft.mp)           # Number of PFTs for multiple panels
n.pft.bp      = length(pft.bp)           # Number of PFTs for box plots by season
n.dbh         = length(dbh.use)          # Number of DBH classes
n.dbh.mp      = length(dbh.mp)           # Number of DBHs for multiple panels
n.dbh.bp      = length(dbh.bp)           # Number of DBHs for box plots
n.season      = length(season.use)       # Number of seasons
n.season.mp   = n.season-1               # Number of seasons for multiple panels
n.year        = length(year.use)         # Number of years
n.drought     = length(drought.use)      # Number of drought lengths
y.sel         = year.use %in% year.avg   # This will select the years for the output
a.pft         = npft+1                   # Index for all PFTs
a.dbh         = ndbh+1                   # Index for all DBH classes
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     Create season suffix.                                                                #
#------------------------------------------------------------------------------------------#
season.suffix = paste(sprintf("%2.2i",sequence(n.season)),tolower(season.key),sep="-")
pft.suffix    = paste("pft",sprintf("%2.2i",pft.use),sep="")
dbh.suffix    = paste("dbh",tolower(dbh.key),sep="-")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      This is the R object that has the simulation information.                           #
#------------------------------------------------------------------------------------------#
rdata.siminfo = file.path( rdata.path
                         , paste( comp.prefix,stext.default,drain.default,"SimInfo.RData"
                                , sep = "_"
                                )#end paste
                         )#end file.path
#------------------------------------------------------------------------------------------#



#----- Retrieve the data. -----------------------------------------------------------------#
cat (" + Loading simulation information from ",basename(rdata.siminfo),"...","\n")
load(rdata.siminfo)
n.simul = simul$n.simul
n.dims  = simul$n.dims
#------------------------------------------------------------------------------------------#


#----- Biomass rate that will cause . -----------------------------------------------------#
crit.change = 100. * log(1.0 - crit.loss) / n.year
#------------------------------------------------------------------------------------------#


#----- Suffix for return period. ----------------------------------------------------------#
optim.suffix = paste( ifelse(nls.optim,ifelse(robust,"nlrob","nls"),"optim")
                    , ifelse(skew.optim & ! nls.optim,"skewn","gauss")
                    , sep = "_"
                    )#end paste
#------------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#
#      Check the loops.                                                                    #
#------------------------------------------------------------------------------------------#
if (is.logical(use.global) && all(use.global)){
   loop.global = sequence(max(1,simul$global$n.level  ))
}else{
   loop.global = use.global
}#end if
n.iphen      = length  (panel$iphen$key   )
n.drain      = length  (scenario$drain$key)
n.stext      = length  (scenario$stext$key)
loop.iphen   = sequence(n.iphen)
loop.stext   = sequence(n.stext)
loop.drought = sequence(n.drought)
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Find the best set up for plotting all seasons and all PFTs in the same plot.        #
#------------------------------------------------------------------------------------------#
lo.season   = pretty.box(n=n.season.mp        )
lo.pft      = pretty.box(n=n.pft.mp           )
lo.dbh      = pretty.box(n=n.dbh.mp           )
lo.iphen    = pretty.box(n=n.iphen            )
if (n.panel == 0 | n.panel == 1){
   lo.panel = pretty.box(n=simul$panel$n.level)
}else if (n.panel == 2){
   lo.panel = pretty.box(n=simul$dim[simul$dim.type == "panel"])
}#end if
#==========================================================================================#
#==========================================================================================#




#------------------------------------------------------------------------------------------#
#     Loop over all global dimensions.                                                     #
#------------------------------------------------------------------------------------------#
for (g in loop.global){


   #---------------------------------------------------------------------------------------#
   #     Select the runs that belong to this global dimension.                             #
   #---------------------------------------------------------------------------------------#
   if (simul$global$n.level == 0){
      g.sel         = rep(TRUE,times=n.simul)
      global.suffix = "comp"
      global.desc   = ""
   }else{
      g.sel         = simul$global$index == g
      global.suffix = simul$global$level[g]
      global.desc   = simul$global$title[g]
   }#end if
   #----- Define useful aliases. ----------------------------------------------------------#
   gindex       = as.data.frame(simul$index[g.sel,])
   g.type       = simul$dim.type == "global"
   n.gsel       = sum(g.sel)
   loop.gsel    = which(g.sel)
   gsel.name    = simul$name[g.sel,]
   real.name    = realisation$key
   drought.name = drought.key
   #---------------------------------------------------------------------------------------#



   #----- Total number of files. ----------------------------------------------------------#
   n.total  = n.gsel * n.real
   RJ.mat   = arrayInd(ind=sequence(n.total),.dim=c(n.real,n.gsel))
   #---------------------------------------------------------------------------------------#





   #----- File name for this global scenario. ---------------------------------------------#
   return.global = file.path( rdata.path
                            , paste("return_sim_",simul$global$level[g],".RData",sep="")
                            )#end file.path
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find out whether to read the return period again.                                #
   #---------------------------------------------------------------------------------------#
   if (file.exists(return.global)){
      #----- File is there, re-open it. ---------------------------------------------------#
      cat("   - Retrieving return period from file: ",basename(return.global),"...","\n")
      load(return.global)
      #------------------------------------------------------------------------------------#
   }else{

      #----- File is not there, start over. -----------------------------------------------#
      cat("   - Starting return period...","\n")
      dft     = list()#end list
      rj.last = 0
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Loop over the data that must be read.                                             #
   #---------------------------------------------------------------------------------------#
   if (rj.last < n.total){
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #     Loop over all files and get the full time series.                              #
      #------------------------------------------------------------------------------------#
      cat ("   - Reading the output from simulations...","\n")
      loop.rj = seq(from=rj.last+1,to=n.total,by=1)
      for (rj in loop.rj){
         #----- Get current indices. ------------------------------------------------------#
         r       = RJ.mat[rj,1]
         j       = RJ.mat[rj,2]
         rj.last = rj
         #---------------------------------------------------------------------------------#




         #------ Grab simulation. ---------------------------------------------------------#
         sim.name  = gsel.name[j,r]
         sim.iphen =                      simul$index[g.sel,"iphen"][j]
         sim.drain = scenario$drain$value[simul$index[g.sel,"drain"][j]]
         sim.stext = scenario$stext$value[simul$index[g.sel,"stext"][j]]
         sim.iata  =                      simul$index[g.sel,"iata" ][j]


         #---------------------------------------------------------------------------------#
         #     Load the data set.                                                          #
         #---------------------------------------------------------------------------------#
         rdata.simul = file.path( here,sim.name,"rdata_month"
                                , paste(sim.name,".RData",sep="")
                                )#end file.path

         cat  ("     * Load data from file ",paste("(",rj,"/",n.total,")",sep="")
                                            ,basename(rdata.simul),"...","\n")
         dummy = load (rdata.simul)
         emean = datum$emean
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the drought length.                                                    #
         #---------------------------------------------------------------------------------#
         dry          = emean$nmon.wdef > 0
         ndry         = length(dry)
         didx         = cumsum(dry & c(FALSE,! dry[-ndry]))
         dry.length   = tapply(X=emean$nmon.wdef    [dry],INDEX=didx[dry],FUN=max )
         dry.begin    = tapply(X=sequence(ndry)     [dry],INDEX=didx[dry],FUN=min )
         dry.end      = tapply(X=sequence(ndry)     [dry],INDEX=didx[dry],FUN=max )
         dry.mwd      = tapply(X=emean$water.deficit[dry],INDEX=didx[dry],FUN=max )
         dry.smpot    = tapply(X=emean$smpot        [dry],INDEX=didx[dry],FUN=max )
         dry.agb.min  = tapply(X=emean$agb          [dry],INDEX=didx[dry],FUN=max )
         dry.agb.mean = tapply(X=emean$agb          [dry],INDEX=didx[dry],FUN=mean)
         dry.ddmort   = tapply(X=emean$agb.ncbmort  [dry],INDEX=didx[dry],FUN=mean)
         dry.dimort   = tapply(X=emean$agb.dimort   [dry],INDEX=didx[dry],FUN=mean)
         dry.agb.1st  = emean$agb[dry.begin]
         dry.agb.last = emean$agb[dry.end  ]
         dry.change   = 1200. * log(dry.agb.last/dry.agb.1st) / (dry.end-dry.begin+1)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find annual rainfall and statistics.                                        #
         #---------------------------------------------------------------------------------#
         yft.rain  = tapply(X=emean$rain, INDEX=datum$year,FUN=sum)
         yft.stat  = sn.stats(yft.rain)
         rain.loc  = yft.stat[1]
         rain.sca  = yft.stat[2]
         rain.shp  = yft.stat[3]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find simulation averages.                                                   #
         #---------------------------------------------------------------------------------#
         sim.agb.1st  = emean$agb[1]
         sim.agb.last = emean$agb[ndry]
         sim.agb.mean = mean(emean$agb)
         sim.agb.min  = min (emean$agb)
         sim.agb.max  = max (emean$agb)
         sim.et       = mean(emean$wflxca*yr.day)
         #---------------------------------------------------------------------------------#




         #------ Bind all data from this simulation to a list element. --------------------#
         dft[[rj]] = list( idx          = rj            + 0 * dry.length
                         , realisation  = r - 1         + 0 * dry.length
                         , iata         = sim.iata      + 0 * dry.length
                         , iphen        = sim.iphen     + 0 * dry.length
                         , drain        = sim.drain     + 0 * dry.length
                         , stext        = sim.stext     + 0 * dry.length
                         , length       = dry.length    + 0 * dry.length
                         , mwd          = dry.mwd       + 0 * dry.length
                         , smpot        = dry.smpot     + 0 * dry.length
                         , agb.min      = dry.agb.min   + 0 * dry.length
                         , agb.mean     = dry.agb.mean  + 0 * dry.length
                         , agb.1st      = dry.agb.1st   + 0 * dry.length
                         , agb.last     = dry.agb.last  + 0 * dry.length
                         , change       = dry.change    + 0 * dry.length
                         , ddmort       = dry.ddmort    + 0 * dry.length
                         , dimort       = dry.dimort    + 0 * dry.length
                         , rain.loc     = rain.loc      + 0 * dry.length
                         , rain.sca     = rain.sca      + 0 * dry.length
                         , rain.shp     = rain.shp      + 0 * dry.length
                         , sim.agb.1st  = sim.agb.1st   + 0 * dry.length
                         , sim.agb.last = sim.agb.last  + 0 * dry.length
                         , sim.agb.mean = sim.agb.mean  + 0 * dry.length
                         , sim.agb.min  = sim.agb.min   + 0 * dry.length
                         , sim.agb.max  = sim.agb.max   + 0 * dry.length
                         , sim.et       = sim.et        + 0 * dry.length
                         )#end list
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Save it if it is time to save.                                               #
         #---------------------------------------------------------------------------------#
         if (rj == n.total){
            #----- Final save, we convert it to data frame. -------------------------------#
            dft = data.frame( apply(X=sapply(X=dft,FUN=c),MARGIN=1,FUN=unlist)
                            , row.names        = NULL
                            , stringsAsFactors = FALSE
                            )#end data.frame
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Map the evapotranspiration of the observed case (r+000, real-00) to all  #
            # cases with the same soil texture and site.                                   #
            #------------------------------------------------------------------------------#
            cat("     * Map the default evapotranspiration onto the default scenario..."
               ,"\n")
            z.sel       = ( dft$drain       == which(scenario$drain$key == "r+000")
                          & dft$realisation == 0 )
            zft         = dft[z.sel,]
            iuse        = apply(X=dft,MARGIN=1,FUN=idx.zrain,default=zft)
            dft$sim.et0 = zft$sim.et[iuse]
            #------------------------------------------------------------------------------#


            #----- Save data to an R object file. -----------------------------------------#
            cat("   - Save return period to: ",basename(return.global),"...","\n")
            save(list=c("rj.last","dft"),file=return.global)
            #------------------------------------------------------------------------------#
         }else if ((rj %% save.every) == 0){
            #----- Save data to an R object file. -----------------------------------------#
            cat("   - Save return period to: ",basename(return.global),"...","\n")
            save(list=c("rj.last","dft"),file=return.global)
            cat("   - Quitting...","\n")
            q("no")
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (rj in loop.rj)
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   }#end if (file.exists(return.global))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find return period for multiple-year droughts.                                   #
   #---------------------------------------------------------------------------------------#
   cat("     * Find drought return periods...","\n")
   my.drought = list()
   pret.lim   = list()
   change.lim = list()
   summ.key   = c("y0","y0.se","y0.p","a","a.se","a.p","b","b.se","b.p"
                 ,"x.crit","x.crit.se","r2")
   n.summ     = length(summ.key)
   summ.tab   = array( data     = NA
                     , dim      = c(n.stext,n.iphen,n.summ,n.drought)
                     , dimnames = list( scenario$stext$key
                                      , panel$iphen$key
                                      , summ.key
                                      , drought.key
                                      )# end dimnames
                     )#end array
   #---------------------------------------------------------------------------------------#


   #----- Loop over drought lengths. ------------------------------------------------------#
   for (d in loop.drought){
      cat("       # >= ",drought.desc[d],"...","\n",sep="")
      my.drought[[d]] = list()
      pret.lim  [[d]] = NULL
      change.lim[[d]] = NULL


      #----- Scale the return period and stress threshold by the number of years. ---------#
      event      = floor ( dft$length / (12 * drought.use[d]) )
      dry.count  = tapply( X = event    , INDEX = dft$idx, FUN = sum   )
      dry.stext  = tapply( X = dft$stext, INDEX = dft$idx, FUN = median)
      dry.iphen  = tapply( X = dft$iphen, INDEX = dft$idx, FUN = median)
      dry.change = ( 100. 
                   * log( tapply( X = dft$sim.agb.last, INDEX = dft$idx, FUN = median)
                        / tapply( X = dft$sim.agb.1st , INDEX = dft$idx, FUN = median) )
                   / n.year )
      n.simul    = length(dry.count)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    In addition, calculate the expected return period for those that did not        #
      # experience any drought.                                                            #
      #------------------------------------------------------------------------------------#
      dry.location = tapply( X = dft$rain.loc, INDEX = dft$idx, FUN = median)
      dry.scale    = tapply( X = dft$rain.sca, INDEX = dft$idx, FUN = median)
      dry.shape    = tapply( X = dft$rain.shp, INDEX = dft$idx, FUN = median)
      dry.et0      = tapply( X = dft$sim.et0 , INDEX = dft$idx, FUN = median)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Randomly sample rainfall to estimate the probabilities.  We then aggregate   #
      # by the number of years and the sample size.                                        #
      #------------------------------------------------------------------------------------#
      if (R.major <= 2){
         random = mapply( FUN      = rsn
                        , location = dry.location
                        , scale    = dry.scale
                        , shape    = dry.shape
                        , MoreArgs = list(n=drought.use[d]*n.sample)
                        )#end mapply
      }else{
         random = mapply( FUN      = rsn
                        , xi       = dry.location
                        , omega    = dry.scale
                        , alpha    = dry.shape
                        , MoreArgs = list(n=drought.use[d]*n.sample)
                        )#end mapply
      }#end if (R.major <= 2)
      random = array(data=random,dim=c(n.sample,drought.use[d],n.simul))
      random = aperm(a = random,perm=c(3,1,2))
      random = apply(X = random, MARGIN=c(1,2),FUN=sum)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Estimate the probability based on the number of successes.  The return-period #
      # is the inverse of the probability, scaled by the number of years so its units are  #
      # always in years.                                                                   #
      #------------------------------------------------------------------------------------#
      thresh.dry    = rep(dry.et0 * drought.use[d],times=n.sample) + 0 * random
      pdr.now       = rowSums(random <= thresh.dry) / n.sample
      zero          = pdr.now == 0.
      pret.now      = drought.use[d] / pdr.now
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     For the simulations that actually experienced the droughts, replace them by    #
      # the actual number of events.                                                       #
      #------------------------------------------------------------------------------------#
      hit           = dry.count > 0
      pret.now[hit] = n.year / dry.count[hit]
      pch.now       = ifelse(dry.count == 0,4,2)
      #------------------------------------------------------------------------------------#



      #------ Update range of return periods (for plots). ---------------------------------#
      pret.lim  [[d]] = range(pret.now,finite=TRUE)
      change.lim[[d]] = range(dry.change,finite=TRUE)
      #------------------------------------------------------------------------------------#



      #---------------------------------------------------------------------------------------#
      #     Loop over all phenologies.                                                        #
      #---------------------------------------------------------------------------------------#
      for (p in loop.iphen){
         #------------------------------------------------------------------------------------#
         #     Select the runs that belong to this phenology.                                 #
         #------------------------------------------------------------------------------------#
         p.sel                 = dry.iphen == p
         my.drought[[d]][[p]] = list()
         cat("         ~ Phenology: ",p," - ",panel$iphen$desc[p],"\n",sep="")
         #------------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Find one fit for each soil texture.                                       #
         #---------------------------------------------------------------------------------#
         for (s in loop.stext){
            #------------------------------------------------------------------------------#
            s.now = scenario$stext$value[s]
            s.sel = dry.stext == s.now
            sel   = p.sel & s.sel & is.finite(pret.now)
            #------------------------------------------------------------------------------#
            cat("           = Soil texture: ",s," - ",scenario$stext$desc[s],"\n",sep="")




            #------------------------------------------------------------------------------#
            #      We define the danger zone by fitting a hyperbolic curve.                #
            #------------------------------------------------------------------------------#
            data.now  = data.frame( change = dry.change[sel]
                                  , pret   = pret.now  [sel]
                                  , pch    = pch.now   [sel]
                                  )#end data.frame
            pret.span = range(log(pret.now[sel]),na.rm=TRUE)
            data.pred = data.frame( pret = exp(seq( from       = pret.span[1]
                                                  , to         = pret.span[2]
                                                  , length.out = n.pred     ) ) )
            #------------------------------------------------------------------------------#



            #----- Use nls to fit the data. -----------------------------------------------#
            if (nls.optim){
               ans = change.return.nls  ( datum   = data.now
                                        , y.crit  = crit.change
                                        , verbose = verbose.optim
                                        , robust  = robust
                                        )#end change.return.optim
            }else{
               ans = change.return.optim( datum   = data.now
                                        , y.crit  = crit.change
                                        , skew    = skew.optim
                                        , verbose = verbose.optim
                                        )#end change.return.optim
            }#end if
            #------------------------------------------------------------------------------#


            #----- Predict the curve. -----------------------------------------------------#
            y0.now           = ans$coefficients[1]
            y0.p             = ans$p.value     [1]
            y0.se            = ans$std.err     [1]
            a.now            = ans$coefficients[2]
            a.p              = ans$p.value     [2]
            a.se             = ans$std.err     [2]
            b.now            = ans$coefficients[3]
            b.p              = ans$p.value     [3]
            b.se             = ans$std.err     [3]
            x.min            = 0.1*min(data.now$pret)
            x.max            = max(data.now$pret)
            x.crit           = ans$x.crit
            x.crit.se        = ans$x.crit.se
            data.pred$change = y0.now + a.now / data.pred$pret^b.now
            data.pred$prob   = drought.use[d] / data.pred$pret
            #------------------------------------------------------------------------------#



            #----- Save predicted data to an array. ---------------------------------------#
            summ.tab[s,p,,d] = c( y0.now, y0.se     , y0.p
                                , a.now , a.se      , a.p
                                , b.now , b.se      , b.p
                                , x.crit, x.crit.se , ans$r.square 
                                )#end c
            #------------------------------------------------------------------------------#




            #----- Save data set. ---------------------------------------------------------#
            my.drought[[d]][[p]][[s]] = list( summ     = ans
                                            , data     = data.now
                                            , pred     = data.pred
                                            )#end list
            #------------------------------------------------------------------------------#
         }#end for (s in loop.stext)
         #---------------------------------------------------------------------------------#
      }#end for (p in loop.iphen)
      #------------------------------------------------------------------------------------#
   }#end for (d in loop.drought)
   #---------------------------------------------------------------------------------------#










   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #     Plot biomass loss as a function of drought length.                                #
   #---------------------------------------------------------------------------------------#
   cat(" + Plotting biomass loss as a function of drought length...","\n",sep="")
   for (p in loop.iphen){
      cat("   - Phenology: ",p," - ",panel$iphen$desc[p],"\n",sep="")
      #----- Retrieve return data. --------------------------------------------------------#
      sel = dft$iphen == p
      rp  = dft[sel,]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Split the drought length into bins.                                            #
      #------------------------------------------------------------------------------------#
      cat("     * Find distribution of biomass change by drought length...","\n")
      brks.dlen   = c(0,seq(from=6,to=quantile(c(rp$length),prob=0.975),by=12),Inf)
      nbrks       = length(brks.dlen)
      labels.dlen = paste (brks.dlen[-nbrks],brks.dlen[-1],sep="-")
      cut.dlen    = cut   ( x      = rp$length
                          , breaks = brks.dlen
                          , right  = FALSE
                          , labels = labels.dlen
                          )#end cut.dlen
      idx.dlen    = match (cut.dlen,levels(cut.dlen))
      #------------------------------------------------------------------------------------#



      #----- Split the change in biomass into drought length bins. ------------------------#
      ok              = ! ( is.na(rp$change) | is.na(rp$stext) | is.na(cut.dlen) )
      myfacts         = list(rp$stext[ok],cut.dlen[ok])
      change.dl       = split(x=rp$change[ok],f=myfacts)
      change.dl.mean  = sapply(X=change.dl,FUN=mean)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find the break points for dry season lengths, leaving enough room for all     #
      # soil textures.                                                                     #
      #------------------------------------------------------------------------------------#
      x.grid  = seq(from=0,by=n.stext+1,length.out=nbrks)
      x.at    = ( rep( x = seq(from=1,by=1        ,length.out=n.stext), times = nbrks-1)
                + rep( x = seq(from=0,by=n.stext+1,length.out=nbrks-1), each  = n.stext)
                )#end x.at
      x.label = parse(text=c(brks.dlen[-nbrks],"infinity"))
      xlimit  = pretty.xylim(u=range(x.grid   ),fracexp=0.0,is.log=FALSE)
      ylimit  = pretty.xylim(u=range(rp$change),fracexp=0.0,is.log=FALSE)
      #------------------------------------------------------------------------------------#



      #------ Set plot annotation. --------------------------------------------------------#
      letitre = paste(panel$iphen$desc[p])
      lex     = desc.unit(desc="Drought length",unit=untab$month   )
      ley     = desc.unit(desc="AGB change"    ,unit=untab$pcagboyr)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Create the box plots with with change in biomass as a function of dry          #
      # season length.                                                                     #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open file or display. -----------------------------------------------------#
         fichier = paste("bp_change_drylen_",panel$iphen$key[p],".",outform[o],sep="")
         fichier = file.path(outroot,fichier)
         if (outform[o] == "x11"){
            X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=wsize$width,height=wsize$height,pointsize=ptsz
                      ,paper=wsize$paper)
         }else if(outform[o] == "pdf"){
            pdf(file=fichier,onefile=FALSE,width=wsize$width,height=wsize$height
               ,pointsize=ptsz,paper=wsize$paper)
         }#end if
         #---------------------------------------------------------------------------------#



         #------ Split the window into two. -----------------------------------------------#
         par(par.user)
         layout(mat=rbind(2,1),heights=c(5,1))
         #---------------------------------------------------------------------------------#




         #------ First plot, the legend. --------------------------------------------------#
         par(mar=c(0.1,4.6,0.1,4.6))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x      = "center"
               , inset  = 0.0
               , legend = paste(scenario$stext$desc," (",scenario$stext$label,")"
                               ,sep="")
               , fill   = scenario$stext$colour
               , border = foreground
               , title  = expression(bold("Soil texture classes"))
               , ncol   = min(3,pretty.box(n.stext)$ncol)
               , xpd    = TRUE
               , cex    = 0.9*cex.ptsz
               )#end legend
         #---------------------------------------------------------------------------------#



         #----- Second plot, the box plots. -----------------------------------------------#
         par(mar=c(5.1,4.6,4.1,2.1))
         plot.new()
         plot.window(xlim=xlimit,ylim=ylimit)
         axis   (side=1,at=x.grid,labels=x.label)
         axis   (side=2,las=1)
         abline (h=axTicks(2),v=x.grid,col=grid.colour,lty="dotted")
         abline (h=0,col=foreground,lty="solid")
         abline (h=crit.change,col=red.mg,lwd=2,lty="dotdash")
         boxplot( x      = change.dl
                , col    = rep(scenario$stext$colour,times=nbrks-1)
                , border = foreground
                , axes   = FALSE
                , add    = TRUE
                , at     = x.at
                )#end boxplot
         points (x=x.at,y=change.dl.mean,pch=21,col=foreground,bg=indigo.mg
                ,cex=1.3,lwd=2)
         title(main=letitre,xlab=lex,ylab=ley)
         box()
         #---------------------------------------------------------------------------------#



         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for (o in 1:nout)
      #------------------------------------------------------------------------------------#
   }#end for (p in loop.iphen)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#









   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #     Plot the biomass loss as a function of return period.                             #
   #---------------------------------------------------------------------------------------#
   cat(" + Plotting biomass loss as a function of return period...","\n",sep="")
   if (robust && nls.optim){
      form.fit = expression(paste("Robust Fit: ",delta[A*G*B] == delta[0] + a / tau[D*1]^b))
   }else{
      form.fit = expression(paste("MLE Fit: ",delta[A*G*B] == delta[0] + a / tau[D*1]^b))
   }#end if

   if (pch.return.type){
      leg.fit = c( expression(paste("MLE Fit: "
                                   ,delta[A*G*B] == delta[0] + a / tau[D*1]^b))
                 , expression(paste("Predicted Critical value: ("
                                   ,tau[c],";",delta[c],")"))
                 , expression(" ")
                 , expression(paste("Event based (",N[D*1] > 0,")"))
                 , expression(paste("Probability based (",N[D*1] == 0,")"))
                 , expression(" ")
                 )#end c
   }else{
      leg.fit = c( expression(paste("MLE Fit: "
                                   ,delta[A*G*B] == delta[0] + a / tau[D*1]^b))
                 , expression(paste("Predicted Critical value: ("
                                   ,tau[c],";",delta[c],")"))
                 , expression(" ")
                 )#end c
   }#end if

   for (d in loop.drought){
      cat("   - ",drought.desc[d],"...","\n",sep="")


      #------------------------------------------------------------------------------------#
      #      Define the range for all drought return periods.                              #
      #------------------------------------------------------------------------------------#
      xlimit = pret.lim  [[d]]
      ylimit = change.lim[[d]]
      if (is.null(std.pret.limit)){
         xlimit  = pretty.xylim(u=xlimit,fracexp=0.0,is.log=TRUE )
      }else{
         xlimit  = c(max(std.pret.limit[1],drought.use[d]),std.pret.limit[2])
      }#end if
      ylimit     = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
      #------------------------------------------------------------------------------------#





      #------ Set plot annotation. --------------------------------------------------------#
      letitre = drought.desc[d]
      lex     = desc.unit(desc="Drought Return period",unit=untab$yr      )
      ley     = desc.unit(desc="Biomass change"       ,unit=untab$pcagboyr)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Create the simulation loss as a function of drought return periods.            #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open file or display. -----------------------------------------------------#
         fichier = paste("return_loss_",drought.key[d],"_",optim.suffix
                        ,".",outform[o],sep="")
         fichier = file.path(outroot,fichier)
         if (outform[o] == "x11"){
            X11(width=size$width,height=size$height,pointsize=ptsz)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=size$width*depth,height=size$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=size$width,height=size$height,pointsize=ptsz
                      ,paper=size$paper)
         }else if(outform[o] == "pdf"){
            pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
               ,pointsize=ptsz,paper=size$paper)
         }#end if
         #---------------------------------------------------------------------------------#



         #------ Split the window into two. -----------------------------------------------#
         par(oma=c(0,2,3,0))
         par(par.user)
         layout(mat     = rbind(rep(lo.iphen$mat.off2,each=2)
                               ,c(rep(1,times=lo.iphen$ncol)
                                 ,rep(2,times=lo.iphen$ncol)
                                 )#end c
                               )#end rbind
               ,heights = c(rep(5/lo.iphen$nrow,times=lo.iphen$nrow),1)
               )#end layout
         #---------------------------------------------------------------------------------#




         #------ First plot, the legend. --------------------------------------------------#
         par(mar=c(0.1,4.1,0.1,0.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         if (pch.return.type){
             legend( x      = "center"
                   , inset  = 0.0
                   , legend = paste(scenario$stext$desc," (",scenario$stext$label,")"
                                   ,sep="")
                   , fill   = scenario$stext$colour
                   , border = foreground
                   , title  = expression(bold("Soil texture classes"))
                   , ncol   = min(2,pretty.box(n.stext)$ncol)
                   , xpd    = TRUE
                   , cex    = 0.8*cex.ptsz
                   )#end legend
         }else{
             legend( x      = "center"
                   , inset  = 0.0
                   , legend = paste(scenario$stext$desc," (",scenario$stext$label,")"
                                   ,sep="")
                   , col    = scenario$stext$colour
                   , pch    = scenario$stext$pch
                   , title  = expression(bold("Soil texture classes"))
                   , ncol   = min(2,pretty.box(n.stext)$ncol)
                   , xpd    = TRUE
                   , cex    = 0.8*cex.ptsz
                   )#end legend
         }#end if
         #---------------------------------------------------------------------------------#




         #------ Second plot, the legend. -------------------------------------------------#
         par(mar=c(0.1,0.1,0.1,2.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x      = "center"
               , inset  = 0.0
               , legend = leg.fit
               , col    = foreground
               , pch    = c(-1,-1,-1,2,4,-1)
               , lwd    = 2
               , pt.lwd = 2
               , lty    = c(1,4,0,0,0,0)
               , title  = expression(bold("Description"))
               , ncol   = ifelse(pch.return.type,2,1)
               , xpd    = TRUE
               , cex    = 0.8*cex.ptsz
               )#end legend
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Loop over the phenologies.                                                 #
         #---------------------------------------------------------------------------------#
         for (p in loop.iphen){
            lacle = paste(panel$iphen$desc[p])

            #----- Find out where the box goes, and set up axes and margins. --------------#
            left    = lo.iphen$left  [p]
            right   = lo.iphen$right [p]
            top     = lo.iphen$top   [p]
            bottom  = lo.iphen$bottom[p]
            mar.now = lo.iphen$mar   [p,] # c(2 + 2 * bottom,1 + 2 * left,0 + 2 * top,0 + 2 * right) + 0.1
            #------------------------------------------------------------------------------#



            #----- Plot the biomass change as a function of drought return period. --------#
            par(mar=mar.now)
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log="x")
            if (bottom) axis(side=1)
            if (left  ) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline (h=0,col=foreground,lty="solid")
            abline (h=crit.change,col=indigo.fg,lwd=2,lty="dotdash")
            for (s in loop.stext){
               dry   = my.drought[[d]][[p]][[s]]
               pred  = dry$pred
               datum = dry$data
               s.col = scenario$stext$colour[s]
               if (pch.return.type){
                  pch.now = datum$pch
               }else{
                  pch.now = scenario$stext$pch   [s]
               }#end if
               abline (v=summ.tab[s,p,"x.crit",d] ,col=s.col,lwd=1      ,lty="dotdash")
               lines  (x=pred$pret ,y=pred$change ,col=s.col,lwd=2      ,lty="solid"  )
               points (x=datum$pret,y=datum$change,col=s.col,pch=pch.now,cex=0.8      )
            }#end for
            title  (main=lacle,cex.main=0.9*cex.main)
            box()
            #------------------------------------------------------------------------------#
         }#end for (d in loop.drought)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot global title.                                                          #
         #---------------------------------------------------------------------------------#
         gtitle( main      = letitre
               , xlab      = lex
               , ylab      = ley
               , off.xlab  = 1/6
               , line.xlab = 4.1
               , line.ylab = 2.6
               , cex.main  = 1.1*cex.ptsz
               )#end gtitle
         #---------------------------------------------------------------------------------#



         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for (o in 1:nout)
      #------------------------------------------------------------------------------------#
   }#end for (d in loop.drought)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


   #------ Re-format table for easy display. ----------------------------------------------#
   summ.tab = signif(summ.tab[,,,1],4)
   idx      = arrayInd( ind       = seq_along(summ.tab)
                      , .dim      = dim(summ.tab)
                      , .dimnames = dimnames(summ.tab)
                      )#end arrayInd
   summ.tab = split(x=summ.tab,f=idx[,1])
   summ.tab = lapply( X   = summ.tab
                    , FUN = matrix
                    , nrow = n.iphen
                    , ncol = n.summ
                    , dimnames = list(panel$iphen$key,summ.key)
                    )#end lapply
   names(summ.tab) = scenario$stext$key
   #---------------------------------------------------------------------------------------#
}#end for (g in loop.global)
#==========================================================================================#
#==========================================================================================#
