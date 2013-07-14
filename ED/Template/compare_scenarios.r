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
here          = getwd()                     # Current directory
srcdir        = "/n/home00/mlongo/util/Rsc" # Script directory
stext.default = "stext16"                   # Default soil texture
drain.default = "r+000"                     # Default rainfall
ibackground   = 0                           # Target background colour:
                                            #    (to adjust foreground colours accordingly)
                                            # 0 -- White
                                            # 1 -- Pitch black
                                            # 2 -- Dark grey
bg.default    = paste("ibg",sprintf("%2.2i",ibackground),sep="")
outroot       = file.path(here
                         ,paste("scencomp",stext.default,drain.default,bg.default,sep="_")
                         )#end file.path
comp.prefix   = "stext"
save.every    = 100                         # Save the partial processing every (save.every)
                                            #    steps
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
retrieve.siminfo = FALSE                            # Retrieve previously loaded sim. info
retrieve.global  = TRUE                             # Retrieve previously loaded data
rdata.path       = file.path(here,"RData_scenario") # Path with R object.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     The following variables control whether to plot the different stuff.  For all of     #
# them, the values are:                                                                    #
#  TRUE  -- plot all of them                                                               #
#  NA    -- plot only the first   (debugging only).                                        #
#  FALSE -- skip them altogether                                                           #
#------------------------------------------------------------------------------------------#
plot.panel       = c(FALSE,TRUE,NA)[2]
plot.tseries     = c(FALSE,TRUE,NA)[1]
plot.szpft       = c(FALSE,TRUE,NA)[2]
plot.barplot     = c(FALSE,TRUE,NA)[2]
plot.xyzvars     = c(FALSE,TRUE,NA)[1]
plot.scencomp    = c(FALSE,TRUE,NA)[2]
plot.panelbox    = c(FALSE,TRUE,NA)[2]
plot.panelxyz    = c(FALSE,TRUE,NA)[2]
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Select the function you would like to use to aggregate years for box plots.  Options #
# are: mean, median, sum (or similar), or in case you want to use all data points (boxplot #
# only, this is not a possibility for others), use unlist or c.                            #
# Important: use the name of the function in quotes.                                       #
#------------------------------------------------------------------------------------------#
tseries.aggr = "median"
boxplot.aggr = "mean"
barplot.aggr = "mean"
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     The following variable allows you to set how to change the limits between runs:      #
# 0 -- Each panel will get the best scale for that particular plot.                        #
# 1 -- All plots are with fixed limits for all panels and scenarios                        #
# 2 -- Plots are with fixed limits for all panels within each fixed scenario               #
# 3 -- Plots are with fixed limits for all panels of the same site for each fixed scenario #
#------------------------------------------------------------------------------------------#
ylim.fix.all = 3
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
                     , colour  = c("#FED164","#FDA531","#FF5308","#D90B00","#6E0500")
                     , pch     = c(12L,13L,5L,6L,8L)
                     , alabel  = c("Soil texture")
                     )#end list
#----- Realisation variables. -------------------------------------------------------------#
key.realisation = sprintf("%2.2i",seq(from=0,to=15,by=1))
realisation     = list( key     = paste("real",key.realisation,sep="-")
                      , desc    = paste("Realisation",key.realisation,sep=" ")
                      , pattern = "real-ZZ"
                      )#end list
#------------------------------------------------------------------------------------------#





#------ Miscellaneous settings. -----------------------------------------------------------#
yeara          = 1952         # First year we will include
yeare          = 1972         # First year to use in the averaged output
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
outform        = c("pdf")              # Formats for output file.  Supported formats are:
                                       #   - "X11" - for printing on screen
                                       #   - "eps" - for postscript printing
                                       #   - "png" - for PNG printing
                                       #   - "pdf" - for PDF printing

byeold         = TRUE                  # Remove old files of the given format?

depth          = 96                    # PNG resolution, in pixels per inch
paper          = "letter"              # Paper size, to define the plot shape
wpaper         = "long"                # Paper size for wide plots.
ptsz           = 16                    # Font size.
lwidth         = 2.5                   # Line width
plotgrid       = TRUE                  # Should I plot the grid in the background? 

legwhere       = "topleft"             # Where should I place the legend?
inset          = 0.01                  # Inset between legend and edge of plot region.
fracexp        = 0.40                  # Expansion factor for y axis (to fit legend)
n.colourbar    = 32                    # Number of colours for the colour bars
n.whitebar     =  1                    # Number of levels around zero to be set to white
notch          = FALSE                 # Add notches to the box plots.
mtext.xoff.im  = -9.00                 # Offset for the x label
mtext.xoff     = -7.50                 # Offset for the x label
mtext.xoff.e   = -4.50                 # Offset for the x label
mtext.xoff.xyz = -3.90                 # Offset for the x label
mtext.yoff     = -0.50                 # Offset for the y label
mtext.yoff.im  = -1.50                 # Offset for the y label
mtext.xadj     =  0.50                 # Offset for the x label
mtext.yadj     =  0.65                 # Offset for the y label
barplot.lwd    =  1.50                 # Line width for the bar plots
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



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#

#----- Set some scaling for objects that may need to shrink for large font sizes. ---------#
cex.ptsz = min(15/ptsz,1.0)
#------------------------------------------------------------------------------------------#


#----- Set some dimensions associated with the simulations. -------------------------------#
n.global      = length(global         )
n.panel       = length(panel          )
n.scenario    = length(scenario       )
n.allscen     = n.scenario + 1
n.realisation = length(realisation$key)
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
   cat (" -- REALISATION =",n.realisation,"\n")
   stop("Correct the number of dimensions or re-write the script...")
}#end if
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size  = plotsize(proje=FALSE,paper=paper)
wsize = plotsize(proje=FALSE,paper=wpaper)
#wsize = modifyList(x=size,val=list(width=size$width*4/3,ratio=size$ratio*4/3))
#------------------------------------------------------------------------------------------#


#----- Make sure that the base directory exists. ------------------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#


#----- Make sure that the RData directory exists. -----------------------------------------#
if (! file.exists(rdata.path)) dir.create(rdata.path)
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
rdata.siminfo = file.path(rdata.path,paste(comp.prefix,"SimInfo.RData",sep="_"))
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#     Here we read or reload the simulations.                                              #
#------------------------------------------------------------------------------------------#
if (retrieve.siminfo && file.exists(rdata.siminfo)){
   #----- Retrieve the data. --------------------------------------------------------------#
   cat (" + Loading simulation information from ",basename(rdata.siminfo),"...","\n")
   load(rdata.siminfo)
   n.simul = simul$n.simul
   n.dims  = simul$n.dims
   #---------------------------------------------------------------------------------------#
}else{
   #---------------------------------------------------------------------------------------#
   #     Load the simulations again.                                                       #
   #---------------------------------------------------------------------------------------#
   cat (" + Reloading the simulation output...","\n")



   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #     In this block we concatenate all dimensions in a way that we can easily modify    #
   # the dimensions later.                                                                 #
   #---------------------------------------------------------------------------------------#
   cat ("   - Finding the parameter space dimensions for these simulations...","\n")
      #----- Combine all the dimensions we are going to explore. --------------------------#
      rdims      = c("global","panel","scenario")
      loop.rdims = seq.len(rdims)
      simul.a    = NULL
      simul.b    = NULL
      dim.type   = NULL
      pattern    = NULL
      default    = NULL
      alabel     = NULL
      for (r in loop.rdims){
         this     = get(rdims[r])
         this.a   = lapply(X=this  ,FUN= data.frame,stringsAsFactors=FALSE)
         this.b   = sapply(X=this.a,FUN=rbind                             )
         this.e   = t(apply(X = sapply(X=this,FUN=c), MARGIN=1,FUN=unlist))
         simul.a  = c    (simul.a,this.a)
         simul.b  = cbind(simul.b,this.b)
         dim.type = c    (dim.type,rep(rdims[r],times=length(this)))

         #---------------------------------------------------------------------------------#
         #    These variables should remain scalars.                                       #
         #---------------------------------------------------------------------------------#
         pattern  = unlist(c(pattern, this.e[,"pattern"]))
         default  = unlist(c(default, this.e[,"default"]))
         alabel   = unlist(c(alabel , this.e[,"alabel" ]))
         #---------------------------------------------------------------------------------#
         rm(this,this.a,this.b,this.e)
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Rename the elements of the 3 simple vectors.                                   #
      #------------------------------------------------------------------------------------#
      names(pattern) = names(default) = names(alabel) = names(simul.a)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find all dimensions and dimension names, then we transform the internal        #
      # vectors in arrays, so it is easier to track                                        #
      #------------------------------------------------------------------------------------#
      dim.simul    = sapply(X=simul.a,FUN=nrow)
      dnames.simul = simul.b    ["key",]
      simul        = apply(X=simul.b,FUN=expand.grid,MARGIN=1,stringsAsFactors=FALSE)
      #------------------------------------------------------------------------------------#



      #----- Append the dimensions to the list. -------------------------------------------#
      simul$dim      = dim.simul
      simul$dimnames = dnames.simul
      simul$dim.type = dim.type
      simul$pattern  = pattern
      simul$default  = default
      simul$alabel   = alabel
      #------------------------------------------------------------------------------------#




      #----- Map the indices onto the dimension names. ------------------------------------#
      simul$index = mapply(FUN= match,x= simul$key,table= simul$dimnames)
      #------------------------------------------------------------------------------------#




      #----- Map the default indices onto the dimension names. ----------------------------#
      simul$def.idx = mapply(FUN= match,x= simul$default,table= simul$dimnames)
      #------------------------------------------------------------------------------------#




      #----- Map the indices onto the dimension names. ------------------------------------#
      simul$name = rep(name.template,times=nrow(simul$index))
      for (l in 1:length(simul$key)){
         simul$name = mapply( FUN         = gsub
                            , replacement = simul$key[[l]]
                            , x           = simul$name
                            , MoreArgs    = list(pattern=simul$pattern[l])
                            )#end mapply
         dimnames(simul$name) = NULL
      }#end for
      #------------------------------------------------------------------------------------#



      #----- Find the global number of runs and dimensions. -------------------------------#
      simul$n.simul = length(simul$name)
      simul$n.dims  = length(simul$dim)
      n.simul       = simul$n.simul
      n.dims        = simul$n.dims
      sim.width     = nchar(n.simul)
      sim.label     = paste( "sim"
                           , sprintf( paste("%",sim.width,".",sim.width,"i",sep="")
                                    , sequence(n.simul))
                           , sep="-"
                           )#end paste


      #------------------------------------------------------------------------------------#
      #     Check whether there are realisations or not.                                   #
      #------------------------------------------------------------------------------------#
      if (n.realisation == 0){
         #---------------------------------------------------------------------------------#
         #     No realisation, still make a dummy matrix.                                  #
         #---------------------------------------------------------------------------------#
         simul$n.real         = 1
         n.real               = simul$n.real
         simul$name           = matrix(simul$name,nrow=n.simul,ncol=n.real)
         dimnames(simul$name) = list(sim.label,"real-00")
         #---------------------------------------------------------------------------------#
      }else{
         #----- Find the combination of all realisations. ---------------------------------#
         simul$n.real = n.realisation
         n.real       = simul$n.real
         simul$name = rep   (simul$name     ,times = n.real )
         replace.by = rep   (realisation$key,each  = n.simul)
         simul$name = matrix( data = mapply( FUN         = gsub
                                           , replacement = replace.by
                                           , x           = simul$name
                                           , MoreArgs    = list(pattern=realisation$pattern)
                                           )#end mapply
                            , nrow = n.simul
                            , ncol = n.real
                            )#end matrix
         dimnames(simul$name) = list(sim.label,realisation$key)
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Create all possible combinations of global and panel variables.               #
      #------------------------------------------------------------------------------------#
      #----- Global variables. ------------------------------------------------------------#
      i.global              = which(simul$dim.type == "global")
      simul$global          = list()
      if (n.global == 1){
         simul$global$index = unlist(simul$key [,i.global])
         simul$global$title = unlist(simul$desc[,i.global])
      }else{
         simul$global$index = apply( X        = simul$key [,i.global]
                                   , MARGIN   = 1
                                   , FUN      = paste
                                   , collapse = "-"
                                   )#end apply
         simul$global$title = apply( X        = simul$desc[,i.global]
                                   , MARGIN   = 1
                                   , FUN      = paste
                                   , collapse = " - "
                                   )#end apply
      }#end if
      simul$global$level      = unique(simul$global$index)
      simul$global$title      = unique(simul$global$title)
      simul$global$index      = match(simul$global$index,simul$global$level)
      simul$global$n.level    = length(simul$global$level)
      #----- Panel variables. -------------------------------------------------------------#
      i.panel               = which(simul$dim.type == "panel")
      simul$panel           = list()
      if (n.panel == 0){
         simul$panel$index  = character(0)
         simul$panel$title  = character(0)
      }else if (n.panel == 1){
         simul$panel$index  = c(simul$key [,i.panel])
         simul$panel$title  = c(simul$desc[,i.panel])
      }else{
         simul$panel$index  = apply( X        = simul$key [,i.panel]
                                   , MARGIN   = 1
                                   , FUN      = paste
                                   , collapse = "-"
                                   )#end apply
         simul$panel$title  = apply( X        = simul$desc[,i.panel]
                                   , MARGIN   = 1
                                   , FUN      = paste
                                   , collapse = " - "
                                   )#end apply
      }#end if
      simul$panel$level         = unique(simul$panel$index)
      simul$panel$title         = unique(simul$panel$title)
      simul$panel$index         = match (simul$panel$index,simul$panel$level)
      simul$panel$n.level       = length(simul$panel$level)
      #====================================================================================#
      #====================================================================================#







      #====================================================================================#
      #====================================================================================#
      #    The scenario variables are done in a slightly different way because when two    #
      # dimensions are given, we must fix one and let the other vary.                      #
      #------------------------------------------------------------------------------------#
      i.scenario   = which(simul$dim.type == "scenario")
      simul$scenario = list()
      if (n.scenario == 1){
         simul$scenario$level  = character(0)
         simul$scenario$title  = character(0)
         simul$scenario$fixcol = integer(0)
         simul$scenario$idxcol = integer(0)
         simul$scenario$index  = list(sequence(n.simul))
         simul$scenario$idxoff = simul$scenario$index
      }else{
         my.keys               = simul$key [,i.scenario]
         my.desc               = simul$desc[,i.scenario]
         un.level              = lapply((X = my.keys),FUN=unique)
         un.title              = lapply((X = my.desc),FUN=unique)
         simul$scenario$level  = unlist(un.level)
         simul$scenario$title  = unlist(un.title)
         simul$scenario$fixcol = rep(i.scenario,times=sapply(X=un.level,FUN=length))
         simul$scenario$idxcol = match(simul$scenario$fixcol,unique(simul$scenario$fixcol))
         simul$scenario$index  = mapply(match,my.keys,un.level)
         #----- Correct the indices so they are in sequence. ------------------------------#
         off                   = c( 0, apply( X       = simul$scenario$index
                                            , MARGIN  = 2
                                            , FUN     = max)[-n.scenario]
                                  )#end c
         names(off)            = colnames(simul$scenario$index)
         simul$scenario$idxoff = ( mapply("+",data.frame(simul$scenario$index),off)
                                 + 0*simul$scenario$index )
      }#end if
      simul$scenario$n.level   = length(simul$scenario$level)
      #------------------------------------------------------------------------------------#


      #----- Delete the temporary arrays. -------------------------------------------------#
      rm(simul.a,simul.b,pattern,default,alabel,dim.simul,dnames.simul)
      #------------------------------------------------------------------------------------#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


   #------ Save the header. ---------------------------------------------------------------#
   cat (" + Saving simulation information from ",basename(rdata.siminfo),"...","\n")
   save(simul,file=rdata.siminfo)
   #---------------------------------------------------------------------------------------#
}#end if
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      Find the best set up for plotting all seasons and all PFTs in the same plot.        #
#------------------------------------------------------------------------------------------#
lo.season   = pretty.box(n=n.season.mp        )
lo.pft      = pretty.box(n=n.pft.mp           )
lo.dbh      = pretty.box(n=n.dbh.mp           )

if (n.panel == 0 | n.panel == 1){
   lo.panel = pretty.box(n=simul$panel$n.level)
}else if (n.panel == 2){
   lo.panel = pretty.box(n=simul$dim[simul$dim.type == "panel"])
}#end if
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Create all directories.                                                             #
#------------------------------------------------------------------------------------------#
outpath = list()
for (o in 1:nout){
   this.out = outform[o]

   onow = list()
   if (! this.out %in% "x11"){
      #------ Main directory for this output format. --------------------------------------#
      onow$main = file.path(outroot,this.out)
      if (! file.exists(onow$main)) dir.create(onow$main)
      #------------------------------------------------------------------------------------#



      #------ Loop over global variables. -------------------------------------------------#
      onow$global = list()
      if (is.logical(use.global) && all(use.global)){
         loop.global    = sequence(max(1,simul$global$n.level  ))
      }else{
         loop.global    = use.global
      }#end if
      loop.panel     = sequence(max(1,simul$panel$n.level   ))
      loop.scenario  = which(simul$scenario$level %in% simul$default)
      loop.allscen   = c(loop.scenario,0)
      for (g in loop.global){
         gnow      = list()

         #---- Main directory for this global combination. --------------------------------#
         gnow$main = file.path(onow$main,simul$global$level[g])
         if (! file.exists(gnow$main)) dir.create(gnow$main)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Create one directory for each scenario, or dump everything in the main     #
         # path if it is a scenario for all.                                               #
         #---------------------------------------------------------------------------------#
         if (n.scenario > 1){
            #------------------------------------------------------------------------------#
            #      Create the directory where the 2-D scenario comparison will be placed.  #
            #------------------------------------------------------------------------------#
            root.scenpanel.season    = file.path(gnow$main,"scenpanel_season")
            gnow$scenpanel.year      = file.path(gnow$main,"scenpanel_year"  )
            root.scenpanelpft.season = file.path(gnow$main,"scenpanelpft_season")
            root.scenpanelpft.year   = file.path(gnow$main,"scenpanelpft_year"  )
            root.scenpaneldbh.season = file.path(gnow$main,"scenpaneldbh_season")
            root.scenpaneldbh.year   = file.path(gnow$main,"scenpaneldbh_year"  )
            if (! file.exists(root.scenpanel.season   )){
               dir.create(root.scenpanel.season   )
            }#end if
            if (! file.exists(root.scenpanelpft.season)){
               dir.create(root.scenpanelpft.season)
            }#end if
            if (! file.exists(root.scenpanelpft.year  )){
               dir.create(root.scenpanelpft.year  )
            }#end if
            if (! file.exists(root.scenpaneldbh.season)){
               dir.create(root.scenpaneldbh.season)
            }#end if
            if (! file.exists(root.scenpaneldbh.year  )){
               dir.create(root.scenpaneldbh.year  )
            }#end if
            if (! file.exists(gnow$scenpanel.year     )){
               dir.create(gnow$scenpanel.year     )
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Create one directory by season.                                         #
            #------------------------------------------------------------------------------#
            gnow$scenpanel.season  = file.path(root.scenpanel.season,season.suffix)
            for (e in season.mp){
               #----- Create one directory for each season. -------------------------------#
               if (! file.exists(gnow$scenpanel.season[e])){
                  dir.create(gnow$scenpanel.season[e])
               }#end if
               #---------------------------------------------------------------------------#

            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Create one directory by PFT or DBH.                                     #
            #------------------------------------------------------------------------------#
            gnow$scenpanelpft.year = file.path(root.scenpanelpft.year,pft.suffix)
            gnow$scenpaneldbh.year = file.path(root.scenpaneldbh.year,dbh.suffix)
            for (f in 1:(n.pft-1)){
               #----- Create one directory for each season. -------------------------------#
               if (! file.exists(gnow$scenpanelpft.year[f])){
                  dir.create(gnow$scenpanelpft.year[f])
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            for (d in 1:n.dbh){
               #----- Create one directory for each season. -------------------------------#
               if (! file.exists(gnow$scenpaneldbh.year[d])){
                  dir.create(gnow$scenpaneldbh.year[d])
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Create one directory by season, then by PFT or DBH.                     #
            #------------------------------------------------------------------------------#
            #----- Combine all PFT/DBHs. and seasons. -------------------------------------#
            pft.season.suffix         = expand.grid(pft.suffix,season.suffix)
            dbh.season.suffix         = expand.grid(dbh.suffix,season.suffix)
            #----- Create the "branches". -------------------------------------------------#
            branch.scenpanelpft.season = apply( X        = pft.season.suffix
                                              , MARGIN   = 1
                                              , FUN      = paste
                                              , collapse ="/"
                                              )#end apply
            branch.scenpaneldbh.season = apply( X        = dbh.season.suffix
                                              , MARGIN   = 1
                                              , FUN      = paste
                                              , collapse ="/"
                                              )#end apply
            gnow$scenpanelpft.season = matrix( data = file.path(root.scenpanelpft.season
                                                               ,branch.scenpanelpft.season)
                                             , nrow = n.pft
                                             , ncol = n.season
                                             )#end matrix
            gnow$scenpaneldbh.season = matrix( data = file.path(root.scenpaneldbh.season
                                                               ,branch.scenpaneldbh.season)
                                             , nrow = n.dbh
                                             , ncol = n.season
                                             )#end matrix
            for (e in season.mp){
               for (f in 1:(n.pft-1)){
                  #----- Create one directory for each PFT. -------------------------------#
                  if (! file.exists(dirname(gnow$scenpanelpft.season[f,e]))){
                     dir.create(dirname(gnow$scenpanelpft.season[f,e]))
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Create one directory for each season. ----------------------------#
                  if (! file.exists(gnow$scenpanelpft.season[f,e])){
                     dir.create(gnow$scenpanelpft.season[f,e])
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               for (d in 1:n.dbh){
                  #----- Create one directory for each DBH. -------------------------------#
                  if (! file.exists(dirname(gnow$scenpaneldbh.season[d,e]))){
                     dir.create(dirname(gnow$scenpaneldbh.season[d,e]))
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Create one directory for each season. ----------------------------#
                  if (! file.exists(gnow$scenpaneldbh.season[d,e])){
                     dir.create(gnow$scenpaneldbh.season[d,e])
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#







            #------------------------------------------------------------------------------#
            #      Create the directory where the multiple box plot will be placed.        #
            #------------------------------------------------------------------------------#
            root.panelbox.season    = file.path(gnow$main,"panelbox_season")
            gnow$panelbox.year      = file.path(gnow$main,"panelbox_year"  )
            root.panelboxpft.season = file.path(gnow$main,"panelboxpft_season")
            root.panelboxpft.year   = file.path(gnow$main,"panelboxpft_year"  )
            root.panelboxdbh.season = file.path(gnow$main,"panelboxdbh_season")
            root.panelboxdbh.year   = file.path(gnow$main,"panelboxdbh_year"  )
            if (! file.exists(root.panelbox.season   )){
               dir.create(root.panelbox.season   )
            }#end if
            if (! file.exists(root.panelboxpft.season)){
               dir.create(root.panelboxpft.season)
            }#end if
            if (! file.exists(root.panelboxpft.year  )){
               dir.create(root.panelboxpft.year  )
            }#end if
            if (! file.exists(root.panelboxdbh.season)){
               dir.create(root.panelboxdbh.season)
            }#end if
            if (! file.exists(root.panelboxdbh.year  )){
               dir.create(root.panelboxdbh.year  )
            }#end if
            if (! file.exists(gnow$panelbox.year     )){
               dir.create(gnow$panelbox.year     )
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Create one directory by season.                                         #
            #------------------------------------------------------------------------------#
            gnow$panelbox.season  = file.path(root.panelbox.season,season.suffix)
            for (e in season.mp){
               #----- Create one directory for each season. -------------------------------#
               if (! file.exists(gnow$panelbox.season[e])){
                  dir.create(gnow$panelbox.season[e])
               }#end if
               #---------------------------------------------------------------------------#

            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Create one directory by PFT or DBH.                                     #
            #------------------------------------------------------------------------------#
            gnow$panelboxpft.year = file.path(root.panelboxpft.year,pft.suffix)
            gnow$panelboxdbh.year = file.path(root.panelboxdbh.year,dbh.suffix)
            for (f in 1:(n.pft-1)){
               #----- Create one directory for each season. -------------------------------#
               if (! file.exists(gnow$panelboxpft.year[f])){
                  dir.create(gnow$panelboxpft.year[f])
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            for (d in 1:n.dbh){
               #----- Create one directory for each season. -------------------------------#
               if (! file.exists(gnow$panelboxdbh.year[d])){
                  dir.create(gnow$panelboxdbh.year[d])
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Create one directory by season, then by PFT or DBH.                     #
            #------------------------------------------------------------------------------#
            #----- Combine all PFT/DBHs. and seasons. -------------------------------------#
            pft.season.suffix         = expand.grid(pft.suffix,season.suffix)
            dbh.season.suffix         = expand.grid(dbh.suffix,season.suffix)
            #----- Create the "branches". -------------------------------------------------#
            branch.panelboxpft.season = apply( X        = pft.season.suffix
                                              , MARGIN   = 1
                                              , FUN      = paste
                                              , collapse ="/"
                                              )#end apply
            branch.panelboxdbh.season = apply( X        = dbh.season.suffix
                                              , MARGIN   = 1
                                              , FUN      = paste
                                              , collapse ="/"
                                              )#end apply
            gnow$panelboxpft.season = matrix( data = file.path(root.panelboxpft.season
                                                               ,branch.panelboxpft.season)
                                             , nrow = n.pft
                                             , ncol = n.season
                                             )#end matrix
            gnow$panelboxdbh.season = matrix( data = file.path(root.panelboxdbh.season
                                                               ,branch.panelboxdbh.season)
                                             , nrow = n.dbh
                                             , ncol = n.season
                                             )#end matrix
            for (e in season.mp){
               for (f in 1:(n.pft-1)){
                  #----- Create one directory for each PFT. -------------------------------#
                  if (! file.exists(dirname(gnow$panelboxpft.season[f,e]))){
                     dir.create(dirname(gnow$panelboxpft.season[f,e]))
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Create one directory for each season. ----------------------------#
                  if (! file.exists(gnow$panelboxpft.season[f,e])){
                     dir.create(gnow$panelboxpft.season[f,e])
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               for (d in 1:n.dbh){
                  #----- Create one directory for each DBH. -------------------------------#
                  if (! file.exists(dirname(gnow$panelboxdbh.season[d,e]))){
                     dir.create(dirname(gnow$panelboxdbh.season[d,e]))
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Create one directory for each season. ----------------------------#
                  if (! file.exists(gnow$panelboxdbh.season[d,e])){
                     dir.create(gnow$panelboxdbh.season[d,e])
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#


               #---- Create the panels comparing the annual means for xyz variables. ------#
               root.panel.xyz  = file.path(gnow$main,"panel_xyz"    )
               if (! file.exists(root.panel.xyz)) dir.create(root.panel.xyz)
               #---------------------------------------------------------------------------#


               #---- Loop over scenarios. -------------------------------------------------#
               gnow$panel.xyz      = matrix(NA_character_,nrow=n.allscen,ncol=npanel.yvar)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all scenarios and Y variables for parameter space and      #
               # create the paths as needed.                                               #
               #---------------------------------------------------------------------------#
               for (s in loop.allscen){
                  odx.s = match(s,loop.allscen)
                  if (s == 0){
                     root.now = file.path(root.panel.xyz,"allscen")
                  }else{
                     root.now = file.path(root.panel.xyz,simul$scenario$level[s])
                  }#end if
                  if (! file.exists(root.now)) dir.create(root.now)
                  for (n in 1:npanel.yvar){
                     gnow$panel.xyz[odx.s,n] = file.path(root.now,panel.xyz$yvar$vname[n])
                     if (! file.exists(gnow$panel.xyz[odx.s,n])){
                        dir.create(gnow$panel.xyz[odx.s,n])
                     }#end if
                  }#end for
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Loop over all possible scenarios.                                          #
         #---------------------------------------------------------------------------------#
         gnow$scenario = list()
         for (s in loop.scenario){
            snow = list()

            #------------------------------------------------------------------------------#
            #     We only create scenario-specific simulations if 2-D scenarios were run.  #
            #------------------------------------------------------------------------------#
            if (n.scenario == 1){
               snow$main = gnow$main
            }else{
               snow$main = file.path(gnow$main,simul$scenario$level[s])
               if (! file.exists(snow$main      )) dir.create(snow$main      )
            }#end if
            #------------------------------------------------------------------------------#


            #----- Create the paths by type of plot. --------------------------------------#
            snow$ts.season          = file.path(snow$main,"ts_season"         )
            snow$ts.year            = file.path(snow$main,"ts_year"           )
            # snow$tspft.season     = See below
            snow$tspft.year         = file.path(snow$main,"tspft_year"        )
            # snow$tsdbh.season     = See below
            snow$tsdbh.year         = file.path(snow$main,"tsdbh_year"        )
            # snow$tspftdbh         = See below
            snow$box.season         = file.path(snow$main,"box_season"        )
            snow$box.year           = file.path(snow$main,"box_year"          )
            snow$boxpft.season      = file.path(snow$main,"boxpft_season"     )
            snow$boxpft.year        = file.path(snow$main,"boxpft_year"       )
            snow$boxdbh.season      = file.path(snow$main,"boxdbh_season"     )
            snow$boxdbh.year        = file.path(snow$main,"boxdbh_year"       )
            snow$barplot.season     = file.path(snow$main,"barplot_season"    )
            snow$barplot.year       = file.path(snow$main,"barplot_year"      )

            root.tspft.season       = file.path(snow$main,"tspft_season"      )
            root.tsdbh.season       = file.path(snow$main,"tsdbh_season"      )
            root.tspftdbh           = file.path(snow$main,"tspftdbh"          )
            root.boxpftdbh          = file.path(snow$main,"boxpftdbh"         )
            root.xyz.season         = file.path(snow$main,"xyz_season"        )
            root.xyz.pft            = file.path(snow$main,"xyz_pft"           )
            root.xyz.dbh            = file.path(snow$main,"xyz_dbh"           )
            root.xyz.all            = file.path(snow$main,"xyz_all"           )

            if (! file.exists(snow$ts.season         )) dir.create(snow$ts.season         )
            if (! file.exists(snow$ts.year           )) dir.create(snow$ts.year           )
            if (! file.exists(snow$tspft.year        )) dir.create(snow$tspft.year        )
            if (! file.exists(snow$tsdbh.year        )) dir.create(snow$tsdbh.year        )
            if (! file.exists(snow$box.season        )) dir.create(snow$box.season        )
            if (! file.exists(snow$box.year          )) dir.create(snow$box.year          )
            if (! file.exists(snow$boxpft.season     )) dir.create(snow$boxpft.season     )
            if (! file.exists(snow$boxpft.year       )) dir.create(snow$boxpft.year       )
            if (! file.exists(snow$boxdbh.season     )) dir.create(snow$boxdbh.season     )
            if (! file.exists(snow$boxdbh.year       )) dir.create(snow$boxdbh.year       )
            if (! file.exists(snow$barplot.season    )) dir.create(snow$barplot.season    )
            if (! file.exists(snow$barplot.year      )) dir.create(snow$barplot.year      )
            if (! file.exists(root.tspft.season      )) dir.create(root.tspft.season      )
            if (! file.exists(root.tsdbh.season      )) dir.create(root.tsdbh.season      )
            if (! file.exists(root.tspftdbh          )) dir.create(root.tspftdbh          )
            if (! file.exists(root.boxpftdbh         )) dir.create(root.boxpftdbh         )
            if (! file.exists(root.xyz.season        )) dir.create(root.xyz.season        )
            if (! file.exists(root.xyz.pft           )) dir.create(root.xyz.pft           )
            if (! file.exists(root.xyz.dbh           )) dir.create(root.xyz.dbh           )
            if (! file.exists(root.xyz.all           )) dir.create(root.xyz.all           )
            #------------------------------------------------------------------------------#



            #----- Directories that have sub-directories. ---------------------------------#
            snow$tspft.season       = file.path(root.tspft.season   ,pft.suffix   )
            snow$tsdbh.season       = file.path(root.tsdbh.season   ,dbh.suffix   )
            snow$tspftdbh           = file.path(root.tspftdbh       ,dbh.suffix   )
            snow$boxpftdbh          = file.path(root.boxpftdbh      ,season.suffix)
            #------------------------------------------------------------------------------#


            #----- Create the sub-subdirectories by PFT. ----------------------------------#
            for (f in sequence(n.pft-1)){
               if (! file.exists(snow$tspft.season[f])) dir.create(snow$tspft.season[f])
            }#end for
            #------------------------------------------------------------------------------#


            #----- Create the sub-subdirectories by PFT. ----------------------------------#
            for (d in sequence(n.dbh)){
               if (! file.exists(snow$tsdbh.season[d])) dir.create(snow$tsdbh.season[d])
               if (! file.exists(snow$tspftdbh    [d])) dir.create(snow$tspftdbh    [d])
            }#end for
            #------------------------------------------------------------------------------#


            #----- Create the sub-subdirectories by season. -------------------------------#
            for (e in sequence(n.season)){
               if (! file.exists(snow$boxpftdbh[e])) dir.create(snow$boxpftdbh[e])
            }#end for
            #------------------------------------------------------------------------------#



            #----- Generate the names of the sub-sub-directories by variable. -------------#
            snow$xyz.season   = file.path(root.xyz.season    ,scen.xyz$yvar$vname)
            snow$xyz.pft      = file.path(root.xyz.pft       ,scen.xyz$yvar$vname)
            snow$xyz.dbh      = file.path(root.xyz.dbh       ,scen.xyz$yvar$vname)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all parameter space variables (Y) and create the paths as     #
            # needed.                                                                      #
            #------------------------------------------------------------------------------#
            for (n in 1:nscen.yvar){
               if (! file.exists(snow$xyz.season[n])) dir.create(snow$xyz.season[n])
               if (! file.exists(snow$xyz.pft   [n])) dir.create(snow$xyz.pft   [n])
               if (! file.exists(snow$xyz.dbh   [n])) dir.create(snow$xyz.dbh   [n])
            }#end for
            #------------------------------------------------------------------------------#



            #----- Initialise the name for the XYZ plots by PFT, DBH, and season. ---------#
            root.xyz.all.yvar = file.path(root.xyz.all       ,scen.xyz$yvar$vname)
            snow$xyz.all      = array(NA_character_,dim=c(n.dbh,n.season,nscen.yvar))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all parameter space variables (Y), DBH classes, and season,   #
            # and create the paths as needed.                                              #
            #------------------------------------------------------------------------------#
            for (n in 1:nscen.yvar){
               if (! file.exists(root.xyz.all.yvar[n])) dir.create(root.xyz.all.yvar[n])
               for (e in 1:n.season){
                  for (d in 1:n.dbh){
                      path.now = paste("dbh_",dbh.key[d],"-season_",season.suffix[e],sep="")
                      snow$xyz.all[d,e,n] = file.path(root.xyz.all.yvar[n],path.now)
                      if (! file.exists(snow$xyz.all[d,e,n])){
                         dir.create(snow$xyz.all[d,e,n])
                      }#end if
                      #--------------------------------------------------------------------#
                  }#end for (d in 1:n.dbh)
                  #------------------------------------------------------------------------#
               }#end for (e in 1:n.season.mp)
               #---------------------------------------------------------------------------#
            }#end for (n in 1:nscen.yvar)
            #------------------------------------------------------------------------------#


           gnow$scenario[[s]] = snow
         }#end (s in sequence(max(1,simul$scenario$n.level))
         #---------------------------------------------------------------------------------#

         onow$global[[g]] = gnow
      }#end for (g in 1:simul$global$n.level)
      #------------------------------------------------------------------------------------#
   }#end if (! this.out %in% "x11")
   #---------------------------------------------------------------------------------------#

   outpath[[o]] = onow
}#end for (o in 1:nout)
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#    For simplicity, detach the results from the structure, and free memory from simul.    #
#------------------------------------------------------------------------------------------#
eft       = simul$eft
simul$eft = NULL
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#      Time to plot.  Here we must loop over all scenarios, nested within the panels,      #
# which are nested within the global.                                                      #
#------------------------------------------------------------------------------------------#
cat(" + Processing data for:","\n",sep="")
if (is.logical(use.global) && all(use.global)){
   loop.global = sequence(max(1,simul$global$n.level  ))
}else{
   loop.global = use.global
}#end if
if (is.na(plot.panel)){
   loop.panel     = sequence(max(1,simul$panel$n.level))[1]
}else if (plot.panel){
   loop.panel     = sequence(max(1,simul$panel$n.level))
}else{
   loop.panel     = integer(0)
}#end if
loop.scenario   = which(simul$scenario$level %in% simul$default)
loop.allscen    = c(loop.scenario,0)
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
   g.type    = simul$dim.type == "global"
   n.gsel    = sum(g.sel)
   loop.gsel = which(g.sel)
   gsel.name = simul$name[g.sel,]
   real.name = realisation$key
   #---------------------------------------------------------------------------------------#



   #----- File name for this global scenario. ---------------------------------------------#
   rdata.global = file.path( rdata.path
                           , paste(comp.prefix,"_sim_",simul$global$level[g],".RData"
                                  ,sep="")
                           )#end file.path
   #---------------------------------------------------------------------------------------#



   #----- File name for this global scenario. ---------------------------------------------#
   rdata.status = file.path( rdata.path
                           , paste("status_",comp.prefix,"_sim_",simul$global$level[g]
                                  ,".txt",sep="")
                           )#end file.path
   #---------------------------------------------------------------------------------------#





   #------ Find the total number of simulations to read. ----------------------------------#
   n.total = n.gsel * n.real
   #---------------------------------------------------------------------------------------#



   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #     Retrieve or load the data.                                                        #
   #---------------------------------------------------------------------------------------#
   if (retrieve.global && file.exists(rdata.global)){
      #----- Grab data from previously loaded variable. -----------------------------------#
      cat("   - Retrieving data from file :",basename(rdata.global),"...","\n")
      load(rdata.global)
      eft = get(simul$global$level[g])
      rm(list=simul$global$level[g])
      n.when     = length(eft$when)
      #------------------------------------------------------------------------------------#




      #------ Check which global simulation was last read. --------------------------------#
      if (rj.last == n.total){
         loop.rj = numeric(0)
      }else{
         loop.rj = seq(from=rj.last+1,to=n.total,by=1)
      }#end if
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #     Make some arrays with the right dimensions.                                    #
      #------------------------------------------------------------------------------------#
      empty          = rep  (NA,times=n.when)
      empty.pft      = array(NA,dim=c(n.when,n.pft))
      empty.pftdbh   = array(NA,dim=c(n.when,n.dbh,n.pft))
      ts.array       = array( data     = NA
                            , dim      = c   (       n.gsel,   n.real,  n.year,  n.season)
                            , dimnames = list(gsel.name[,1],real.name,year.key,season.key)
                            )#end array
      tspft.array    = array( data     = NA
                            , dim      = c   (        n.gsel,   n.real,  n.year
                                             ,      n.season,    n.pft)
                            , dimnames = list( gsel.name[,1],real.name,year.key
                                             ,    season.key,  pft.key)
                            )#end array
      tspftdbh.array = array( data     = NA
                            , dim      = c   (       n.gsel,   n.real,  n.year
                                             ,     n.season,    n.dbh,  n.pft)
                            , dimnames = list(gsel.name[,1],real.name,year.key
                                             ,   season.key,  dbh.key,pft.key)
                            )#end array
      #------------------------------------------------------------------------------------#
   }else{
      cat("   - Start loading data from global: ",simul$global$title[g],"...","\n")

      eft = list()

      #------------------------------------------------------------------------------------#
      #     Create a general time stamp that works for all simulations.                    #
      #------------------------------------------------------------------------------------#
      #----- Get all times. ---------------------------------------------------------------#
      eft$when  = c(chron(paste(12,1,yeara-1,sep="/"))
                   ,chron(paste(rep(1:12,times=n.year),1,rep(year.use,each=12),sep="/")))
      eft$when  = eft$when[-length(eft$when)]
      eft$year  = numyears (eft$when)
      eft$month = nummonths(eft$when)
      n.when    = length(eft$when)
      #---- Find the seasons. -------------------------------------------------------------#
      eft$season    = season(eft$when,add.year=TRUE,dec.next=TRUE)
      eft$ss.year   = as.numeric(substring(eft$season,1,4))
      eft$ss.season = as.numeric(substring(eft$season,5,6))
      #---- Find unique identifiers for seasons and years based on season, not months. ----#
      eft$toseason  = unique(eft$season)
      eft$toyear    = unique(eft$ss.year)
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #     Make some arrays with the right dimensions.                                    #
      #------------------------------------------------------------------------------------#
      empty          = rep  (NA,times=n.when)
      empty.pft      = array(NA,dim=c(n.when,n.pft))
      empty.pftdbh   = array(NA,dim=c(n.when,n.dbh,n.pft))
      ts.array       = array( data     = NA
                            , dim      = c   (       n.gsel,   n.real,  n.year,  n.season)
                            , dimnames = list(gsel.name[,1],real.name,year.key,season.key)
                            )#end array
      tspft.array    = array( data     = NA
                            , dim      = c   (        n.gsel,   n.real,  n.year
                                             ,      n.season,    n.pft)
                            , dimnames = list( gsel.name[,1],real.name,year.key
                                             ,    season.key,  pft.key)
                            )#end array
      tspftdbh.array = array( data     = NA
                            , dim      = c   (       n.gsel,   n.real,  n.year
                                             ,     n.season,    n.dbh,  n.pft)
                            , dimnames = list(gsel.name[,1],real.name,year.key
                                             ,   season.key,  dbh.key,pft.key)
                            )#end array
      #------------------------------------------------------------------------------------#






      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #     Initialise all variables.                                                      #
      #------------------------------------------------------------------------------------#
      for (v in 1:nscen.ts){
         #----- Copy variable info. -------------------------------------------------------#
         var.vname  = scen.ts$vname[v]
         var.desc   = scen.ts$desc [v]
         var.quant  = scen.ts$quant[v]
         is.pft     = scen.ts$pft  [v]
         is.dbh     = scen.ts$dbh  [v]
         var.pft    = paste(var.vname,"pft"   ,sep="")
         var.pftdbh = paste(var.vname,"pftdbh",sep="")

         #----- Append the lists to a common name. ----------------------------------------#
         eft[[var.vname]]                      = list()
         eft[[var.vname]]$ts                   = ts.array
         if (is.pft) eft[[var.vname]]$tspft    = tspft.array
         if (is.dbh) eft[[var.vname]]$tspftdbh = tspftdbh.array
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:nscen.ts)
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


      #------ Set the counter to the beginning. -------------------------------------------#
      rj.last = 0
      loop.rj = seq(from=rj.last+1,to=n.total,by=1)
      #------------------------------------------------------------------------------------#

   }#end if
   #---------------------------------------------------------------------------------------#





   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #      Loop over all global variables.                                                  #
   #---------------------------------------------------------------------------------------#
   RJ.mat  = arrayInd(ind=sequence(n.total),.dim=c(n.real,n.gsel))
   for (rj in loop.rj){
      #----- Get current indices. ---------------------------------------------------------#
      r       = RJ.mat[rj,1]
      j       = RJ.mat[rj,2]
      rj.last = rj
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Load the data set.                                                             #
      #------------------------------------------------------------------------------------#
      rdata.simul = file.path(here,gsel.name[j,r],"rdata_month"
                             ,paste(gsel.name[j,r],"RData",sep="."))
      if (file.exists(rdata.simul)){
         cat  ("     * Load data from file ",paste("(",rj,"/",n.total,")",sep="")
                                            ,basename(rdata.simul),"...","\n")
         dummy = try(load(rdata.simul))

         #------ Make sure that the data have been properly loaded, otherwise, quit. ------#
         if ("try-error" %in% is(dummy)){
            cat("       ~ Failed loading data.  Check input!","\n")

            #----- Save status file so we know it is complete. ----------------------------#
            cat ("       ~ Saving report on corruped file to "
                ,basename(rdata.status),"...","\n")
            status = data.frame( rj = rj, total = n.total, failed = basename(rdata.simul))
            write.table( x         = status
                       , file      = rdata.status
                       , append    = FALSE
                       , quote     = FALSE
                       , row.names = FALSE
                       , col.names = TRUE
                       )#end write.table
            cat(" + Quitting...","\n")
            q("no")
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         emean = datum$emean
         szpft = datum$szpft
         #---------------------------------------------------------------------------------#
      }else{
         cat ("       ~ File ",basename(rdata.simul)," doesn't exist!","\n")

         #----- Save status file so we know it is complete. -------------------------------#
         cat ("       ~ Saving report on missing file to ",basename(rdata.status)
             ,"...","\n")
         status = data.frame( rj = rj, total = n.total, failed = basename(rdata.simul))
         write.table( x         = status
                    , file      = rdata.status
                    , append    = FALSE
                    , quote     = FALSE
                    , row.names = FALSE
                    , col.names = TRUE
                    )#end write.table
         cat(" + Quitting...","\n")
         q("no")
      }#end if
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Find the indices for mapping the data from the original data set to the        #
      # combined one.                                                                      #
      #------------------------------------------------------------------------------------#
      idx   = match(eft$when,datum$when)
      w.sel = is.finite(idx)
      idx   = idx[w.sel]
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      These are shorter versions of the season indices.                             #
      #------------------------------------------------------------------------------------#
      ee = sequence(n.season.mp)
      e5 = n.season
      #------------------------------------------------------------------------------------#



      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Copy the time variables to the consolidated list.                              #
      #------------------------------------------------------------------------------------#
      for (v in 1:nscen.ts){
         #----- Copy variable info. -------------------------------------------------------#
         var.vname  = scen.ts$vname [v]
         var.desc   = scen.ts$desc  [v]
         var.f.aggr = get(scen.ts$f.aggr[v])
         var.add    = scen.ts$add   [v]
         var.mult   = scen.ts$mult  [v]
         is.pft     = scen.ts$pft   [v]
         is.dbh     = scen.ts$dbh   [v]
         is.mort    = scen.ts$mort  [v]
         is.recr    = scen.ts$recr  [v]
         cat  ("       > Processing ",var.desc,"...","\n")
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Grab the time series.                                                       #
         #---------------------------------------------------------------------------------#
         var.now        = empty
         if (is.pft){
            var.now[w.sel] = szpft[[var.vname]][idx,a.dbh,a.pft]
         }else{
            var.now[w.sel] = emean[[var.vname]][idx]
         }#end if
         #------ Change units if required by the user. ------------------------------------#
         var.now = var.add + var.now * var.mult
         #------ Find the means/sum by year and by season. --------------------------------#
         eft[[var.vname]]$ts[j,r,,ee] = tapply( X     = var.now
                                              , INDEX = list(eft$ss.year
                                                            ,eft$ss.season)
                                              , FUN   = var.f.aggr
                                              , na.rm = TRUE
                                              )#end tapply
         eft[[var.vname]]$ts[j,r,,e5] = tapply( X     = var.now
                                              , INDEX = eft$ss.year
                                              , FUN   = var.f.aggr
                                              , na.rm = TRUE
                                              )#end tapply
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check what to do depending on whether the variable is a PFT and/or          #
         # DBH.                                                                            #
         #---------------------------------------------------------------------------------#
         if (is.pft){
            var.now         = empty.pft
            var.now[w.sel,] = szpft[[var.vname]][idx,a.dbh,pft.use]
            #------ Change units if required by the user. ---------------------------------#
            var.now = var.add + var.now * var.mult


            #----- Find the means/sum by year and by season. ------------------------------#
            eft[[var.vname]]$tspft[j,r,,ee,] = qapply( X     = var.now
                                                     , INDEX = list( eft$ss.year
                                                                   , eft$ss.season
                                                                   )#end list
                                                     , DIM   = 1
                                                     , FUN   = var.f.aggr
                                                     , na.rm = TRUE
                                                     )#end qapply
            eft[[var.vname]]$tspft[j,r,,e5,] = qapply( X     = var.now
                                                     , INDEX = eft$ss.year
                                                     , DIM   = 1
                                                     , FUN   = var.f.aggr
                                                     , na.rm = TRUE
                                                     )#end qapply
            #------------------------------------------------------------------------------#
         }#end if (var.pft)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check what to do depending on whether the variable is a PFT and/or          #
         # DBH.                                                                            #
         #---------------------------------------------------------------------------------#
         if (is.dbh){
            var.now          = empty.pftdbh
            var.now[w.sel,,] = szpft[[var.vname]][idx,dbh.use,pft.use]
            #------ Change units if required by the user. ---------------------------------#
            var.now = var.add + var.now * var.mult

            #----- Find the means/sum by year and by season. ------------------------------#
            eft[[var.vname]]$tspftdbh[j,r,,ee,,] = qapply( X     = var.now
                                                         , INDEX = list(eft$ss.year
                                                                       ,eft$ss.season
                                                                       )#end list
                                                         , DIM   = 1
                                                         , FUN   = var.f.aggr
                                                         , na.rm = TRUE
                                                         )#end qapply
            eft[[var.vname]]$tspftdbh[j,r,,e5,,] = qapply( X     = var.now
                                                         , INDEX = eft$ss.year
                                                         , DIM   = 1
                                                         , FUN   = var.f.aggr
                                                         , na.rm = TRUE
                                                         )#end qapply
            #------------------------------------------------------------------------------#
         }#end if (var.dbh)
         #---------------------------------------------------------------------------------#




         #----- Transform mortality data. -------------------------------------------------#
         if (is.mort){
            mort.ts          = 100. * (1. - exp(-eft[[var.vname]]$ts[j,r,,]))
            eft[[var.vname]]$ts         [j,r,,  ] = mort.ts
            rm(mort.ts)
            #----- PFT-dependent variable. ------------------------------------------------#
            if (is.pft){
               mort.tspft    = 100. * (1. - exp(-eft[[var.vname]]$tspft[j,r,,,]))
               eft[[var.vname]]$tspft   [j,r,,, ] = mort.tspft
               rm(mort.tspft)
            }#end if
            #------------------------------------------------------------------------------#


            #----- DBH-dependent variable. ------------------------------------------------#
            if (is.dbh){
               mort.tspftdbh = 100. * (1. - exp(-eft[[var.vname]]$tspftdbh[j,r,,,,]))
               eft[[var.vname]]$tspftdbh[j,r,,,,] = mort.tspftdbh
               rm(mort.tspftdbh)
            }#end if
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Transform recruitment data. -----------------------------------------------#
         if (is.recr){
            recr.ts          = 100. * (exp(eft[[var.vname]]$ts[j,r,,]) - 1.0)
            eft[[var.vname]]$ts         [j,r,,  ] = recr.ts
            rm(recr.ts)
            #----- PFT-dependent variable. ------------------------------------------------#
            if (is.pft){
               recr.tspft    = 100. * (exp(eft[[var.vname]]$tspft[j,r,,,]) - 1.0)
               eft[[var.vname]]$tspft   [j,r,,, ] = recr.tspft
               rm(recr.tspft)
            }#end if
            #------------------------------------------------------------------------------#

            #----- DBH-dependent variable. ------------------------------------------------#
            if (is.dbh){
               recr.tspftdbh = 100. * (exp(eft[[var.vname]]$tspftdbh[j,r,,,,]) - 1.0)
               eft[[var.vname]]$tspftdbh[j,r,,,,] = recr.tspftdbh
               rm(recr.tspftdbh)
            }#end if (is.dbh)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:nscen.ts)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





      #------------------------------------------------------------------------------------#
      #     Remove the temporary variables.                                                #
      #------------------------------------------------------------------------------------#
      rm(datum,emean,szpft)
      #------------------------------------------------------------------------------------#



      #----- Save object simul to the R file but keep using eft. --------------------------#
      if ((rj %% save.every) == 0 || rj == n.total){


         #----- Save this variable to the global structure. -------------------------------#
         cat (" + Copying simulation to the structure...","\n")
         dummy = assign(x=simul$global$level[g],value=eft)
         #---------------------------------------------------------------------------------#



         #----- Save partial data frame to RData and free memory. -------------------------#
         cat (" + Saving the simulation output to ",basename(rdata.global),"...","\n")
         save(list=c("rj.last",simul$global$level[g]),file=rdata.global)
         rm  (list=simul$global$level[g])
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      R becomes slow once too many files are loaded.  Quit if not all files      #
         # have been processed, so it refreshes the memory.                                #
         #---------------------------------------------------------------------------------#
         if (rj != n.total){
            cat (" + Quitting...","\n")
            q("no")
         }else{

            #----- Save status file so we know it is complete. ----------------------------#
            cat (" + Saving the simulation status to ",basename(rdata.status),"...","\n")
            status = data.frame( rj = rj, total = n.total, complete = rj == n.total)
            write.table( x         = status
                       , file      = rdata.status
                       , append    = FALSE
                       , quote     = FALSE
                       , row.names = FALSE
                       , col.names = TRUE
                       )#end write.table
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for (rj in loop.rj)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#







   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #     Loop over the panel combinations.                                                 #
   #---------------------------------------------------------------------------------------#
   for (p in loop.panel){
      #------------------------------------------------------------------------------------#
      #     Select the runs that belong to this panel.                                     #
      #------------------------------------------------------------------------------------#
      if (simul$panel$n.level == 0){
         p.sel = rep(TRUE,times=n.gsel)
         cat("     * Single panel","\n",sep="")
      }else{
         p.sel = simul$panel$index[g.sel] == p
         cat("     * Panel: ",p," - ",simul$panel$title[p],"\n",sep="")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find out the site we are processing now, and select all simulations for this   #
      # site.                                                                              #
      #------------------------------------------------------------------------------------#
      iata.now = unique(simul$key$iata[g.sel][p.sel])
      i.sel    = simul$key$iata[g.sel] == iata.now
      #------------------------------------------------------------------------------------#



      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #       Loop over all scenarios.                                                     #
      #------------------------------------------------------------------------------------#
      for (s in loop.scenario){
         #---------------------------------------------------------------------------------#
         #     Select the runs that belong to this panel.                                  #
         #---------------------------------------------------------------------------------#
         if (simul$scenario$n.level == 0){
            s.sel = rep(TRUE,times=n.gsel)
            cat("       > Single scenario","\n",sep="")
         }else{
            s.sel = simul$scenario$idxoff[g.sel,simul$scenario$idxcol[s]] == s
            cat("       > Scenario:",simul$scenario$title[s],"\n",sep="")
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Decide which simulations to include for plot limit estimation.              #
         #---------------------------------------------------------------------------------#
         sel   = p.sel & s.sel
         if (ylim.fix.all == 0){
            l.sel = sel
         }else if (ylim.fix.all == 1){
            l.sel = rep(TRUE,times=n.gsel)
         }else if (ylim.fix.all == 2){
            l.sel = s.sel
         }else if (ylim.fix.all == 3){
            l.sel = i.sel & s.sel
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Build the total selection, the total description, and the total suffix.    #
         #---------------------------------------------------------------------------------#
         if (simul$scenario$n.level == 0){
            if (simul$panel$n.level == 0){
               if (simul$global$n.level == 0){
                  #----- No panels or global.  Only scenarios. ----------------------------#
                  out.desc   = NULL
                  out.suffix = "scencomp"
                  #------------------------------------------------------------------------#
               }else{
                  #----- Scenarios and global. --------------------------------------------#
                  out.desc   = simul$global$title[g]
                  out.suffix = simul$global$level[g]
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#
            }else{
               if (simul$global$n.level == 0){
                  #----- Scenarios and panel. ---------------------------------------------#
                  out.desc   = simul$panel$title[p]
                  out.suffix = simul$panel$level[p]
                  #------------------------------------------------------------------------#
               }else{
                  #----- Scenarios, panel, and global. ------------------------------------#
                  out.desc   = paste(simul$panel$title[p],simul$global$title[g],sep=" - ")
                  out.suffix = paste(simul$panel$level[p],simul$global$level[g],sep="-")
                  #------------------------------------------------------------------------#
               }#end if (simul$global$n.level == 0)
               #---------------------------------------------------------------------------#
            }#end if (simul$panel$n.level == 0)
            #------------------------------------------------------------------------------#
         }else{
            if (simul$panel$n.level == 0){
               if (simul$global$n.level == 0){
                  #----- 2-D scenarios but no panels or global variables. -----------------#
                  out.desc   = simul$scenario$title[s]
                  out.suffix = simul$scenario$level[s]
                  #------------------------------------------------------------------------#
               }else{
                  #----- 2-D scenarios and global variables. ------------------------------#
                  out.desc   = paste(simul$scenario$title[s]
                                    ,simul$global$title  [g]
                                    ,sep=" - ")
                  out.suffix = paste(simul$scenario$level[s]
                                    ,simul$global$level  [g]
                                    ,sep="-")
                  #------------------------------------------------------------------------#
               }#end if (simul$global$n.level == 0)
               #---------------------------------------------------------------------------#
            }else{
               if (simul$global$n.level == 0){
                  #----- 2-D scenarios and panel variables. -------------------------------#
                  out.desc   = paste(simul$scenario$title[s]
                                    ,simul$panel$title   [p]
                                    ,sep=" - ")
                  out.suffix = paste(simul$scenario$level[s]
                                    ,simul$panel$level   [p]
                                    ,sep="-")
                  #------------------------------------------------------------------------#
               }else{
                  #----- 2-D scenarios, panel, and global variables. ----------------------#
                  out.desc   = paste(paste(simul$scenario$title[s]
                                          ,simul$panel$title   [p]
                                          ,sep=" - ")
                                    ,simul$global$title  [g]
                                    ,sep="\n")
                  out.suffix = paste(simul$scenario$level[s]
                                    ,simul$panel$level   [p]
                                    ,simul$global$level  [g]
                                    ,sep="-")
                  #------------------------------------------------------------------------#
               }#end if (simul$global$n.level == 0)
               #---------------------------------------------------------------------------#
            }#end if (simul$panel$n.level == 0)
            #------------------------------------------------------------------------------#
         }#end if (simul$scenario$n.level == 0)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Find the column with the plot information.                                   #
         #---------------------------------------------------------------------------------#
         if (simul$scenario$n.level == 0){
            now = scenario[[1]]
         }else{
            now = scenario[-simul$scenario$idxcol[s]][[1]]
         }#end if
         n.now = length(now$key)
         n.lax = n.now + 2
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the loop for all plots.                                                #
         #---------------------------------------------------------------------------------#
         if (is.na(plot.tseries)){
            loop.tseries = sequence(nscen.ts)[1]
         }else if (plot.tseries){
            loop.tseries = sequence(nscen.ts)
         }else{
            loop.tseries = integer(0)
         }#end if
         if (is.na(plot.szpft)){
            loop.szpft = sequence(nscen.szpft)[1]
         }else if (plot.szpft){
            loop.szpft = sequence(nscen.szpft)
         }else{
            loop.szpft = integer(0)
         }#end if
         if (is.na(plot.barplot)){
            loop.barplot = sequence(nscen.barplot)[1]
         }else if (plot.barplot){
            loop.barplot = sequence(nscen.barplot)
         }else{
            loop.barplot = integer(0)
         }#end if
         if (is.na(plot.xyzvars)){
            loop.x = sequence(nscen.xvar)[1]
            loop.y = sequence(nscen.yvar)[1]
            loop.z = sequence(nscen.zvar)[1]
         }else if (plot.xyzvars){
            loop.x = sequence(nscen.xvar)
            loop.y = sequence(nscen.yvar)
            loop.z = sequence(nscen.zvar)
         }else{
            loop.x = integer(0)
            loop.y = integer(0)
            loop.z = integer(0)
         }#end if
         #---------------------------------------------------------------------------------#


         #.................................................................................#
         #.................................................................................#
         #     Plot the comparison between time series, by season.                         #
         #---------------------------------------------------------------------------------#
         cat  ("         ~ Plotting the time series...","\n")
         for (v in loop.tseries){
            #----- Copy variable info. ----------------------------------------------------#
            var.vname    = scen.ts$vname   [v]
            var.desc     = scen.ts$desc    [v]
            var.unit     = scen.ts$unit    [v]
            is.pft       = scen.ts$pftvar  [v]
            is.dbh       = scen.ts$dbhvar  [v]
            var.plog     = scen.ts$plog    [v]
            var.plog.dbh = scen.ts$plog.dbh[v]
            var.plt      = scen.ts$plt     [v]
            if (var.plog    ){var.xylog      = "y"}else{var.xylog      = ""}
            if (var.plog.dbh){var.xylog.dbh  = "y"}else{var.xylog.dbh  = ""}
            #------------------------------------------------------------------------------#


            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #     Decide whether to plot or not.                                           #
            #------------------------------------------------------------------------------#
            if (var.plt){
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the variable by season.                                         #
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find the X and Y limits for all data (same scale for all).            #
               #---------------------------------------------------------------------------#
               var.now       = eft[[var.vname]]$ts
               every         = apply( X      = var.now
                                    , MARGIN = c(1,3,4)
                                    , FUN    = get(tseries.aggr)
                                    , na.rm  = TRUE
                                    )#end apply
               xlimit        = pretty.xylim( u       = year.use
                                           , fracexp = 0.0
                                           , is.log  = FALSE
                                           )#end pretty.xylim
               ylimit.season = pretty.xylim( u       = every[l.sel,,1:n.season.mp]
                                           , fracexp = 0.0
                                           , is.log   = var.plog
                                           )#end pretty.xylim
               ylimit.year   = pretty.xylim( u       = every[l.sel,,n.season]
                                           , fracexp = 0.0
                                           , is.log  = var.plog
                                           )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #----- Load variable. ------------------------------------------------------#
               this    = every[sel,,]
               ylimit  = ylimit.season
               #---------------------------------------------------------------------------#


               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Seasonal means)","\n",out.desc,sep="")
               lex     = paste("Year")
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$ts.season
                  #------------------------------------------------------------------------#



                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix,"-ts_season."
                                 ,outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(oma = c(0.2,3,4.5,0))
                  layout(mat    = rbind(lo.season$mat.off,rep(1,times=lo.season$ncol))
                        ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
                  par(par.user)
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , col     = now$colour
                         , lwd     = 2.0
                         , pch     = 16
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (e in 1:n.season.mp){
                     #----- Find out where the box goes, and set up axes and margins. -----#
                     left    = (e %% lo.season$ncol) == 1
                     right   = (e %% lo.season$ncol) == 0
                     top     = e <= lo.season$ncol
                     bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(season.desc[e],sep="")
                     #---------------------------------------------------------------------#



                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                     if (bottom) axis(side=1)
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,line=0.5)
                     if (plotgrid) grid(col=grid.colour,lty="solid")

                     #----- Plot the lines. -----------------------------------------------#
                     for (n in 1:n.now){
                        points(x=year.use,y=this[n,,e],type="o",pch=16,col=now$colour[n]
                              ,lwd=2.0)
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.season.mp)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 1/6
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  #par(las=0)
                  #mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  #mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  #mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #     Now we plot the annual data.                                          #
               #---------------------------------------------------------------------------#


               #----- Update limits. ------------------------------------------------------#
               ylimit  = ylimit.year
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means) ","\n",out.desc,sep="")
               lex     = paste("Year")
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath = outpath[[o]]$global[[g]]$scenario[[s]]$ts.year
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix,"-ts_year."
                                 ,outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  layout(mat=rbind(2,1),heights=c(5,1))
                  #------------------------------------------------------------------------#



                  #---- Add the legend. ---------------------------------------------------#
                  par.now = modifyList(x=par.user,val=list(mar=c(0.2,4.0,0.1,2.0)))
                  par(par.now)
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = inset
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , col     = now$colour
                         , lwd     = 2.0
                         , pch     = 16
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=c(5,4,4,2)+0.4)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  axis(side=1)
                  axis(side=2)
                  box()
                  title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.ptsz)
                  if (plotgrid) grid(col=grid.colour,lty="solid")
                  #------------------------------------------------------------------------#



                  #----- Plot the lines. --------------------------------------------------#
                  for (n in 1:n.now){
                     points(x=year.use,y=this[n,,n.season],type="o",pch=16
                           ,col=now$colour[n],lwd=2.0)
                  }#end for
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (var.plt)
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#








            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #     If this is a PFT variable, plot the time series by PFT and season.       #
            #------------------------------------------------------------------------------#
            if (var.plt && is.pft){



               #---------------------------------------------------------------------------#
               #     Grab all data (to build the same scale for all).                      #
               #---------------------------------------------------------------------------#
               var.now       = eft[[var.vname]]$tspft
               every         = apply( X      = var.now
                                    , MARGIN = c(1,3,4,5)
                                    , FUN    = get(tseries.aggr)
                                    , na.rm  = TRUE
                                    )#end apply
               xlimit        = pretty.xylim( u       = year.use
                                           , fracexp = 0.0
                                           , is.log  = FALSE
                                           )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #     Plot the seasonal time series for each PFT.                           #
               #---------------------------------------------------------------------------#
               for (f in 1:(n.pft-1)){



                  #----- Load variable. ---------------------------------------------------#
                  this    = every[sel,,,]
                  #------------------------------------------------------------------------#


                  #----- Y limits for this PFT. -------------------------------------------#
                  ylimit = pretty.xylim( u       = every[l.sel,,-n.season,f]
                                       , fracexp = 0.0
                                       , is.log   = var.plog
                                       )#end pretty.xylim
                  #------------------------------------------------------------------------#



                  #----- Set the title. ---------------------------------------------------#
                  letitre = paste(var.desc," (Seasonal means) - ",pft.desc[f]
                                 ,"\n",out.desc,sep="")
                  lex     = paste("Year")
                  ley     = desc.unit(desc=var.desc,unit=var.unit)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all formats.                                            #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Get the path. -------------------------------------------------#
                     now.outpath = outpath[[o]]$global[[g]]$scenario[[s]]$tspft.season[f]
                     #---------------------------------------------------------------------#


                     #----- Open file or display. -----------------------------------------#
                     fichier = paste(now.outpath,"/",var.vname,"-",pft.suffix[f],"-"
                                    ,out.suffix,"-ts_season.",outform[o],sep="")
                     if (outform[o] == "x11"){
                        X11(width=size$width,height=size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=size$width*depth
                           ,height=size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=size$width,height=size$height
                                  ,pointsize=ptsz,paper=size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                           ,pointsize=ptsz,paper=size$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Split the window into several smaller windows.  Add a bottom    #
                     # row to fit the legend.                                              #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     par.orig = par(no.readonly = TRUE)
                     mar.orig = par.orig$mar
                     par(oma = c(0.2,3,4.5,0))
                     layout(mat    = rbind(lo.season$mat.off,rep(1,times=lo.season$ncol))
                           ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                           )#end layout
                     #---------------------------------------------------------------------#



                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.2,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = if (any(now$parse)){
                                           parse(text=now$legend)
                                        }else{
                                           now$legend
                                        }#end if
                            , col     = now$colour
                            , lwd     = 2.0
                            , pch     = 16
                            , bg      = background
                            , ncol    = min(3,pretty.box(n.now)$ncol)
                            , title   = expression(bold("Simulation"))
                            , cex     = cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Loop over all seasons, and plot the bar plots.                  #
                     #---------------------------------------------------------------------#
                     for (e in 1:n.season.mp){
                        #------------------------------------------------------------------#
                        #    Find out where is this box going, and set up axes and         #
                        # margins.                                                         #
                        #------------------------------------------------------------------#
                        left    = (e %% lo.season$ncol) == 1
                        right   = (e %% lo.season$ncol) == 0
                        top     = e <= lo.season$ncol
                        bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                        #------------------------------------------------------------------#


                        #----- Set up the title for each plot. ----------------------------#
                        lesub = paste(season.desc[e],sep="")
                        #------------------------------------------------------------------#



                        #----- Plot window and grid. --------------------------------------#
                        par(mar=mar.now)
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                        if (bottom) axis(side=1)
                        if (left  ) axis(side=2)
                        box()
                        title(main=lesub,cex.main=cex.ptsz,line=0.5)
                        if (plotgrid) grid(col=grid.colour,lty="solid")

                        #----- Plot the lines. --------------------------------------------#
                        for (n in 1:n.now){
                           points(x=year.use,y=this[n,,e,f],type="o",pch=16
                                 ,col=now$colour[n],lwd=2.0)
                        }#end for
                        #------------------------------------------------------------------#
                     }#end for (e in 1:n.season.mp)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the global title.                                          #
                     #---------------------------------------------------------------------#
                     gtitle( main      = letitre
                           , xlab      = lex
                           , ylab      = ley
                           , off.xlab  = 1/6
                           , line.xlab = 4.1
                           , line.ylab = 2.6
                           , cex.main  = 1.1*cex.ptsz
                           )#end gtitle
                     #par(las=0)
                     #mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                     #mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                     #mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     clean.tmp()
                     #---------------------------------------------------------------------#
                  }#end for (o in 1:nout)
                  #------------------------------------------------------------------------#
               }#end for (f in 1:(n.pft-1))
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the time series by PFT (annual, one panel for each PFT).        #
               #---------------------------------------------------------------------------#
               #----- Set limits. ---------------------------------------------------------#
               ylimit = pretty.xylim( u       = every[l.sel,, n.season,pft.mp]
                                    , fracexp = 0.0
                                    , is.log  = var.plog
                                    )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means)","\n",out.desc,sep="")
               lex     = paste("Year")
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath = outpath[[o]]$global[[g]]$scenario[[s]]$tspft.year
                  #------------------------------------------------------------------------#

                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix,"-tspft_year."
                                 ,outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par.now  = modifyList(x=par.user,val=list(oma = c(0.2,3,4.5,0)))
                  par(par.now)
                  layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                        ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , col     = now$colour
                         , lwd     = 2.0
                         , pch     = 16
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (f in 1:n.pft.mp){

                     #---------------------------------------------------------------------#
                     #     Find out where is this box going, and set up axes and margins.  #
                     #---------------------------------------------------------------------#
                     left    = (f %% lo.pft$ncol) == 1
                     right   = (f %% lo.pft$ncol) == 0
                     top     = f <= lo.pft$ncol
                     bottom  = f > (lo.pft$nrow - 1) * lo.pft$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(pft.desc[pft.mp[f]],sep="")
                     #---------------------------------------------------------------------#



                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                     if (bottom) axis(side=1)
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     if (plotgrid) grid(col=grid.colour,lty="solid")

                     #----- Plot the lines. -----------------------------------------------#
                     for (n in 1:n.now){
                        points(x=year.use,y=this[n,,n.season,pft.mp[f]],type="o",pch=16
                              ,col=now$colour[n],lwd=2.0)
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.pft)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 1/6
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (is.pft && var.plt)
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#








            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #     If this is a DBH variable, plot the time series by DBH and season.       #
            #------------------------------------------------------------------------------#
            if (var.plt && is.dbh){



               #---------------------------------------------------------------------------#
               #     Grab all data (to build the same scale for all).                      #
               #---------------------------------------------------------------------------#
               var.now = apply( X      = eft[[var.vname]]$tspftdbh
                              , MARGIN = c(1,3,4,5,6)
                              , FUN    = get(tseries.aggr)
                              , na.rm  = TRUE
                              )#end apply
               every   = var.now[,,,,n.pft]
               xlimit  = pretty.xylim( u       = year.use
                                     , fracexp = 0.0
                                     , is.log  = FALSE
                                     )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #     Plot the seasonal time series for each PFT.                           #
               #---------------------------------------------------------------------------#
               for (d in 1:n.dbh){



                  #----- Load variable. ---------------------------------------------------#
                  this    = every[sel,,,]
                  #------------------------------------------------------------------------#


                  #----- Y limits for this PFT. -------------------------------------------#
                  ylimit = pretty.xylim( u       = every[l.sel,,-n.season,d]
                                       , fracexp = 0.0
                                       , is.log   = var.plog.dbh
                                       )#end pretty.xylim
                  #------------------------------------------------------------------------#



                  #----- Set the title. ---------------------------------------------------#
                  letitre = paste(var.desc," (Seasonal means) - DBH class: ",dbh.desc[d]
                                 ,"\n",out.desc,sep="")
                  lex     = paste("Year")
                  ley     = desc.unit(desc=var.desc,unit=var.unit)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all formats.                                            #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Get the path. -------------------------------------------------#
                     now.outpath = outpath[[o]]$global[[g]]$scenario[[s]]$tsdbh.season[d]
                     #---------------------------------------------------------------------#


                     #----- Open file or display. -----------------------------------------#
                     fichier = paste(now.outpath,"/",var.vname,"-",dbh.suffix[d],"-"
                                    ,out.suffix,"-ts_season.",outform[o],sep="")
                     if (outform[o] == "x11"){
                        X11(width=size$width,height=size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=size$width*depth
                           ,height=size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=size$width,height=size$height
                                  ,pointsize=ptsz,paper=size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                           ,pointsize=ptsz,paper=size$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Split the window into several smaller windows.  Add a bottom    #
                     # row to fit the legend.                                              #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     par.orig = par(no.readonly = TRUE)
                     mar.orig = par.orig$mar
                     par(oma = c(0.2,3,4.5,0))
                     layout(mat    = rbind(lo.season$mat.off,rep(1,times=lo.season$ncol))
                           ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                           )#end layout
                     #---------------------------------------------------------------------#



                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.2,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = if (any(now$parse)){
                                           parse(text=now$legend)
                                        }else{
                                           now$legend
                                        }#end if
                            , col     = now$colour
                            , lwd     = 2.0
                            , pch     = 16
                            , bg      = background
                            , ncol    = min(3,pretty.box(n.now)$ncol)
                            , title   = expression(bold("Simulation"))
                            , cex     = cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Loop over all seasons, and plot the bar plots.                  #
                     #---------------------------------------------------------------------#
                     for (e in 1:n.season.mp){
                        #------------------------------------------------------------------#
                        #    Find out where is this box going, and set up axes and         #
                        # margins.                                                         #
                        #------------------------------------------------------------------#
                        left    = (e %% lo.season$ncol) == 1
                        right   = (e %% lo.season$ncol) == 0
                        top     = e <= lo.season$ncol
                        bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                        #------------------------------------------------------------------#


                        #----- Set up the title for each plot. ----------------------------#
                        lesub = paste(season.desc[e],sep="")
                        #------------------------------------------------------------------#



                        #----- Plot window and grid. --------------------------------------#
                        par(mar=mar.now)
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n"
                                   ,log=var.xylog.dbh)
                        if (bottom) axis(side=1)
                        if (left  ) axis(side=2)
                        box()
                        title(main=lesub,cex.main=cex.ptsz,line=0.5)
                        if (plotgrid) grid(col=grid.colour,lty="solid")

                        #----- Plot the lines. --------------------------------------------#
                        for (n in 1:n.now){
                           points(x=year.use,y=this[n,,e,d],type="o",pch=16
                                 ,col=now$colour[n],lwd=2.0)
                        }#end for
                        #------------------------------------------------------------------#
                     }#end for (e in 1:n.season.mp)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the global title.                                          #
                     #---------------------------------------------------------------------#
                     gtitle( main      = letitre
                           , xlab      = lex
                           , ylab      = ley
                           , off.xlab  = 1/6
                           , line.xlab = 4.1
                           , line.ylab = 2.6
                           , cex.main  = 1.1*cex.ptsz
                           )#end gtitle
                     # par(las=0)
                     # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                     # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                     # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     clean.tmp()
                     #---------------------------------------------------------------------#
                  }#end for (o in 1:nout)
                  #------------------------------------------------------------------------#
               }#end for (d in 1:n.dbh)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the time series by DBH (annual, one panel for each DBH).        #
               #---------------------------------------------------------------------------#
               #----- Set limits. ---------------------------------------------------------#
               ylimit = pretty.xylim( u       = every[l.sel,, n.season,dbh.mp]
                                    , fracexp = 0.0
                                    , is.log  = var.plog.dbh
                                    )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means)","\n",out.desc,sep="")
               lex     = paste("Year")
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath = outpath[[o]]$global[[g]]$scenario[[s]]$tsdbh.year
                  #------------------------------------------------------------------------#

                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix,"-tsdbh_year."
                                 ,outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par.now  = modifyList(x=par.user,val=list(oma = c(0.2,3,4.5,0)))
                  par(par.now)
                  layout(mat    = rbind(lo.dbh$mat.off,rep(1,times=lo.dbh$ncol))
                        ,height = c(rep(5/lo.dbh$nrow,times=lo.dbh$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , col     = now$colour
                         , lwd     = 2.0
                         , pch     = 16
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (d in 1:n.dbh.mp){

                     #---------------------------------------------------------------------#
                     #     Find out where is this box going, and set up axes and margins.  #
                     #---------------------------------------------------------------------#
                     left    = (d %% lo.dbh$ncol) == 1
                     right   = (d %% lo.dbh$ncol) == 0
                     top     = d <= lo.dbh$ncol
                     bottom  = d > (lo.dbh$nrow - 1) * lo.dbh$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(dbh.mp.desc[d],sep="")
                     #---------------------------------------------------------------------#



                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n"
                                ,log=var.xylog.dbh)
                     if (bottom) axis(side=1)
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     if (plotgrid) grid(col=grid.colour,lty="solid")

                     #----- Plot the lines. -----------------------------------------------#
                     for (n in 1:n.now){
                        points(x=year.use,y=this[n,,n.season,dbh.mp[d]],type="o",pch=16
                              ,col=now$colour[n],lwd=2.0)
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.pft)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 1/6
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (is.pft && var.plt)
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#






            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            if (is.pft && is.dbh && var.plt){
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the time series by DBH and PFT.  Annual, one plot for each DBH, #
               # and within each DBH, one panel for each PFT).                             #
               #---------------------------------------------------------------------------#
               every   = var.now[,,n.season,,]
               #---------------------------------------------------------------------------#


               #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
               #     Loop over DBH classes.                                                #
               #---------------------------------------------------------------------------#
               for (d in 1:n.dbh){
                  #----- Set limits. ------------------------------------------------------#
                  ylimit = pretty.xylim( u       = every[l.sel,,d,pft.mp]
                                       , fracexp = 0.0
                                       , is.log  = var.plog
                                       )#end pretty.xylim
                  #------------------------------------------------------------------------#



                  #----- Set the title. ---------------------------------------------------#
                  letitre = paste(var.desc," (Annual means) - DBH Class: ",dbh.desc[d]
                                 ,"\n",out.desc,sep="")
                  lex     = paste("Year")
                  ley     = desc.unit(desc=var.desc,unit=var.unit)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all formats.                                            #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Get the path. -------------------------------------------------#
                     now.outpath = outpath[[o]]$global[[g]]$scenario[[s]]$tspftdbh
                     #---------------------------------------------------------------------#

                     #----- Open file or display. -----------------------------------------#
                     fichier = paste(now.outpath,"/",var.vname,"-",dbh.suffix[d],"-"
                                    ,out.suffix,"-ts_year.",outform[o],sep="")
                     if (outform[o] == "x11"){
                        X11(width=size$width,height=size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=size$width*depth
                           ,height=size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=size$width,height=size$height
                                  ,pointsize=ptsz,paper=size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                           ,pointsize=ptsz,paper=size$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Split the window into several smaller windows.  Add a bottom    #
                     # row to fit the legend.                                              #
                     #---------------------------------------------------------------------#
                     par.orig = par(no.readonly = TRUE)
                     mar.orig = par.orig$mar
                     par.now  = modifyList(x=par.user,val=list(oma = c(0.2,3,4.5,0)))
                     par(par.now)
                     layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                           ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                           )#end layout
                     #---------------------------------------------------------------------#



                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.2,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = if (any(now$parse)){
                                           parse(text=now$legend)
                                        }else{
                                           now$legend
                                        }#end if
                            , col     = now$colour
                            , lwd     = 2.0
                            , pch     = 16
                            , bg      = background
                            , ncol    = min(3,pretty.box(n.now)$ncol)
                            , title   = expression(bold("Simulation"))
                            , cex     = cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Loop over all seasons, and plot the bar plots.                  #
                     #---------------------------------------------------------------------#
                     for (f in 1:n.pft.mp){

                        #------------------------------------------------------------------#
                        #     Find out where the box is going, set up axes and margins.    #
                        #------------------------------------------------------------------#
                        left    = (f %% lo.pft$ncol) == 1
                        right   = (f %% lo.pft$ncol) == 0
                        top     = f <= lo.pft$ncol
                        bottom  = f > (lo.pft$nrow - 1) * lo.pft$ncol
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                        #------------------------------------------------------------------#


                        #----- Set up the title for each plot. ----------------------------#
                        lesub = paste(pft.mp.desc[f],sep="")
                        #------------------------------------------------------------------#



                        #----- Plot window and grid. --------------------------------------#
                        par(mar=mar.now)
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n",log=var.xylog)
                        if (bottom) axis(side=1)
                        if (left  ) axis(side=2)
                        box()
                        title(main=lesub,cex.main=cex.ptsz,line=0.5)
                        if (plotgrid) grid(col=grid.colour,lty="solid")

                        #----- Plot the lines. --------------------------------------------#
                        for (n in 1:n.now){
                           points(x=year.use,y=every[n,,d,pft.mp[f]],type="o",pch=16
                                 ,col=now$colour[n],lwd=2.0)
                        }#end for
                        #------------------------------------------------------------------#
                     }#end for (f in 1:n.pft.mp)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the global title.                                          #
                     #---------------------------------------------------------------------#
                     gtitle( main      = letitre
                           , xlab      = lex
                           , ylab      = ley
                           , off.xlab  = 1/6
                           , line.xlab = 4.1
                           , line.ylab = 2.6
                           , cex.main  = 1.1*cex.ptsz
                           )#end gtitle
                     # par(las=0)
                     # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                     # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                     # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     clean.tmp()
                     #---------------------------------------------------------------------#
                  }#end for (o in 1:nout)
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               }#end for (d in 1:ndbh)
               #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
            }#end if (is.pft && var.plt)
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
         }#end for (v in 1:nscen.ts){
         #.................................................................................#
         #.................................................................................#






         #.................................................................................#
         #.................................................................................#
         #      Plot the box plots.                                                        #
         #---------------------------------------------------------------------------------#
         cat  ("         ~ Plotting the box plots...","\n")
         for (v in loop.szpft){
            #----- Copy variable info. ----------------------------------------------------#
            var.vname    = scen.szpft$vname   [v]
            var.desc     = scen.szpft$desc    [v]
            var.unit     = scen.szpft$unit    [v]
            is.pft       = scen.szpft$pftvar  [v]
            is.dbh       = scen.szpft$dbhvar  [v]
            var.plog     = scen.szpft$plog    [v]
            var.plog.dbh = scen.szpft$plog.dbh[v]
            var.plt      = scen.szpft$plt     [v]
            if (var.plog    ){var.xylog      = "y"}else{var.xylog      = ""}
            if (var.plog.dbh){var.xylog.dbh  = "y"}else{var.xylog.dbh  = ""}
            #------------------------------------------------------------------------------#






            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            if (var.plt){
               #---------------------------------------------------------------------------#
               #      Plot the box plot by season.                                         #
               #---------------------------------------------------------------------------#
               xlimit = c(0,n.now) + 0.5
               xat    = sequence(n.now)
               xgrid  = seq(from=0.5, to=n.now + 0.5,by=1)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Load all simulations for y axis.                                      #
               #---------------------------------------------------------------------------#
               if (n.real == 1){
                  every = eft[[var.vname]]$ts[,1,y.sel,]
               }else if(boxplot.aggr %in% c("c","unlist")){
                  original = eft[[var.vname]]$ts[,,y.sel,]
                  concat   = apply( X      = expand.grid( dimnames(original)[[2]]
                                                        , dimnames(original)[[3]]
                                                        )#end expand.grid
                                  , MARGIN = 1
                                  , FUN    = paste
                                  , collapse = "+"
                                  )#end apply
                  every    = apply( X      = original
                                  , MARGIN = c(1,4)
                                  , FUN    = unlist
                                  )#end apply
                  every    = aperm(a=every,perm=c(2,1,3))
                  dimnames(every)[[2]] = concat
               }else{
                  every = apply( X      = eft[[var.vname]]$ts[,,y.sel,]
                               , MARGIN = c(1,2,4)
                               , FUN    = get(boxplot.aggr)
                               , na.rm  = TRUE
                               )#end apply
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Get the range for all simulations (fixed y coordinates.               #
               #---------------------------------------------------------------------------#
               ylimit.season = pretty.xylim( u       = every[l.sel,,-n.season]
                                           , fracexp = 0.0
                                           , is.log  = var.plog
                                           )#end pretty.xylimit
               ylimit.year   = pretty.xylim( u       = every[l.sel,, n.season]
                                           , fracexp = 0.0
                                           , is.log  = var.plog
                                           )#end pretty.xylimit
               #---------------------------------------------------------------------------#



               #----- Load variable. ------------------------------------------------------#
               this   = every[sel,,-n.season]
               ylimit = ylimit.season
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Seasonal means)","\n",out.desc,sep="")
               lex     = now$alabel
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the box plot colour.                                           #
               #---------------------------------------------------------------------------#
               bp.colour = now$colour
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$box.season
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-boxplot_season.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wsize$width,height=wsize$height
                               ,pointsize=ptsz,paper=wsize$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=wsize$width,height=wsize$height
                        ,pointsize=ptsz,paper=wsize$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(par.user)
                  par(oma = c(0.2,3,4.5,0))
                  layout(mat    = lo.season$mat)
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (e in 1:n.season.mp){

                     #---------------------------------------------------------------------#
                     #      Find out where is this box going, and set up axes and margins. #
                     #---------------------------------------------------------------------#
                     left    = (e %% lo.season$ncol) == 1
                     right   = (e %% lo.season$ncol) == 0
                     top     = e <= lo.season$ncol
                     bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(season.desc[e],sep="")
                     #---------------------------------------------------------------------#



                     #----- Create the temporary box plot list. ---------------------------#
                     ss.this = this[,,e]
                     ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                     bp.plot = split(x=ss.this,f=list(ai[,1]))
                     #---------------------------------------------------------------------#




                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                     if (bottom) axis(side=1,at=xat,labels=now$label)
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     if (plotgrid){
                        abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                     }#end if

                     #----- Add the box plot, without the x axis. -------------------------#
                     boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                            ,yaxt="n")
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.season.mp)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 0.0
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # mtext(text=lex,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff/3)
                  # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the box plot (annual means).                                    #
               #---------------------------------------------------------------------------#
               xlimit = c(0,n.now) + 0.5
               xat    = sequence(n.now)
               xgrid  = seq(from=0.5, to=n.now + 0.5,by=1)
               #---------------------------------------------------------------------------#


               #----- Load variable and set y axis. ---------------------------------------#
               this   = every[sel,,n.season]
               ylimit = ylimit.year
               #---------------------------------------------------------------------------#




               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means)","\n",out.desc,sep="")
               lex     = now$alabel
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the box plot colour.                                           #
               #---------------------------------------------------------------------------#
               bp.colour = now$colour
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$box.year
                  #------------------------------------------------------------------------#

                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-boxplot_year.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Create the temporary box plot list. ------------------------------#
                  ai      = arrayInd(sequence(length(this)),.dim=dim(this))
                  bp.plot = split(x=this,f=list(ai[,1]))
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(par.user)
                  par(mar=c(5.1,4.4,4.1,2.1))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  axis(side=1,at=xat,labels=now$label)
                  axis(side=2)
                  box()
                  title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.ptsz)
                  if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                  #----- Add the box plot, without the x axis. ----------------------------#
                  boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                         ,yaxt="n")
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (var.plt)
            #..............................................................................#
            #..............................................................................#






            #..............................................................................#
            #..............................................................................#
            if (var.plt && is.pft){
               #---------------------------------------------------------------------------#
               #      Plot the box plot by PFT and by season.                              #
               #---------------------------------------------------------------------------#
               xlimit = c(0,n.pft.bp*n.lax) + 0.5
               xat    = seq(from=1 + 0.5*(n.lax-1), to=n.pft.bp*n.lax, by=n.lax)
               xgrid  = seq(from=0.5, to=(n.pft.bp+1)*n.lax - 0.5*(n.lax-1),by=n.lax)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Load all simulations for y axis.                                      #
               #---------------------------------------------------------------------------#
               if (n.real == 1){
                  every = eft[[var.vname]]$tspft[,1,y.sel,,pft.bp]
               }else if (boxplot.aggr %in% c("c","unlist")){
                  original = eft[[var.vname]]$tspft[,,y.sel,,pft.bp]
                  concat   = apply( X      = expand.grid( dimnames(original)[[2]]
                                                        , dimnames(original)[[3]]
                                                        )#end expand.grid
                                  , MARGIN = 1
                                  , FUN    = paste
                                  , collapse = "+"
                                  )#end apply
                  every    = apply( X      = original
                                  , MARGIN = c(1,4,5)
                                  , FUN    = unlist
                                  )#end apply
                  every    = aperm(a=every,perm=c(2,1,3,4))
                  dimnames(every)[[2]] = concat
               }else{
                  every = apply( X      = eft[[var.vname]]$tspft[,,y.sel,,pft.bp]
                               , MARGIN = c(1,2,4,5)
                               , FUN    = get(boxplot.aggr)
                               , na.rm  = TRUE
                               )#end apply
               }#end if
               ylimit.season = pretty.xylim( u       = every[l.sel,,-n.season,]
                                           , fracexp = 0.0
                                           , is.log  = var.plog
                                           )#end pretty.xylim
               ylimit.year   = pretty.xylim( u       = every[l.sel,, n.season,]
                                           , fracexp = 0.0
                                           , is.log  = var.plog
                                           )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #----- Load variable and set y axis. ---------------------------------------#
               this   = every[sel,,-n.season,]
               ylimit = ylimit.season
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Seasonal means)","\n",out.desc,sep="")
               lex     = paste("Plant functional type")
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the box plot colour.                                           #
               #---------------------------------------------------------------------------#
               bp.colour = rep( x     = c("transparent",now$colour,"transparent")
                              , times = n.pft.bp
                              )#end rep
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$boxpft.season
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-boxpft_season.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wsize$width,height=wsize$height
                               ,pointsize=ptsz,paper=wsize$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=wsize$width,height=wsize$height
                        ,pointsize=ptsz,paper=wsize$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(oma = c(0.2,3,4.5,0))
                  layout(mat    = rbind(lo.season$mat.off,rep(1,times=lo.season$ncol))
                        ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , fill    = now$colour
                         , border  = foreground
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (e in 1:n.season.mp){

                     #---------------------------------------------------------------------#
                     #      Find out where is this box going, and set up axes and margins. #
                     #---------------------------------------------------------------------#
                     left    = (e %% lo.season$ncol) == 1
                     right   = (e %% lo.season$ncol) == 0
                     top     = e <= lo.season$ncol
                     bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(season.desc[e],sep="")
                     #---------------------------------------------------------------------#



                     #----- Create the temporary box plot list. ---------------------------#
                     ss.this = this[,,e,]
                     empty   = array( data     = NA
                                    , dim      = dim(ss.this)[-1]
                                    , dimnames = dimnames(ss.this)[-1]
                                    )#end array
                     ss.this = abind(aaa=empty,ss.this,zzz=empty,along=1)
                     ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                     bp.plot = split(x=ss.this,f=list(ai[,1],ai[,3]))
                     #---------------------------------------------------------------------#




                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                     if (bottom) axis(side=1,at=xat,labels=pft.key[pft.bp])
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")

                     #----- Add the box plot, without the x axis. -------------------------#
                     boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                            ,yaxt="n")
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.season.mp)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 1/6
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the box plot by PFT (annual means).                             #
               #---------------------------------------------------------------------------#
               xlimit = c(0,n.pft.bp*n.lax) + 0.5
               xat    = seq(from=1 + 0.5*(n.lax-1), to=n.pft.bp*n.lax, by=n.lax)
               xgrid  = seq(from=0.5, to=(n.pft.bp+1)*n.lax - 0.5*(n.lax-1),by=n.lax)
               #---------------------------------------------------------------------------#


               #----- Load variable and set y axis. ---------------------------------------#
               this   = every[sel,, n.season,]
               ylimit = ylimit.year
               #---------------------------------------------------------------------------#




               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means)","\n",out.desc,sep="")
               lex     = paste("Plant functional type")
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the box plot colour.                                           #
               #---------------------------------------------------------------------------#
               bp.colour = rep(c("transparent",now$colour,"transparent"),times=n.pft.bp)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$boxpft.year
                  #------------------------------------------------------------------------#

                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-boxpft_year.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Create the temporary box plot list. ------------------------------#
                  empty   = array( data     = NA
                                 , dim      = dim(this)[-1]
                                 , dimnames = dimnames(this)[-1]
                                 )#end array
                  ss.this = abind(aaa=empty,this,zzz=empty,along=1)
                  ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                  bp.plot = split(x=ss.this,f=list(ai[,1],ai[,3]))
                  #------------------------------------------------------------------------#


                  #----- Split into two plots. --------------------------------------------#
                  par(par.user)
                  layout(mat=matrix(c(2,1),nrow=2),heights=c(5,1))
                  #------------------------------------------------------------------------#



                  #----- Add the legend. --------------------------------------------------#
                  par(mar=c(0.2,4.4,0.1,2.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = inset
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , fill    = now$colour
                         , border  = foreground
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=c(3.1,4.4,4.1,2.1))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  axis(side=1,at=xat,labels=pft.key[pft.bp])
                  axis(side=2)
                  box()
                  title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.ptsz)
                  if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                  #----- Add the box plot, without the x axis. ----------------------------#
                  boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                         ,yaxt="n")
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (var.plt)
            #..............................................................................#
            #..............................................................................#






            #..............................................................................#
            #..............................................................................#
            if (var.plt & is.dbh){
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the box plot by DBH and by season.                              #
               #---------------------------------------------------------------------------#
               xlimit = c(0,n.dbh.bp*n.lax) + 0.5
               xat    = seq(from=1 + 0.5*(n.lax-1), to=n.dbh.bp*n.lax, by=n.lax)
               xgrid  = seq(from=0.5, to=(n.dbh.bp+1)*n.lax - 0.5*(n.lax-1), by=n.lax)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Load all simulations for y axis.                                      #
               #---------------------------------------------------------------------------#
               if (n.real == 1){
                  every = eft[[var.vname]]$tspftdbh[,1,y.sel,,dbh.bp,n.pft]
               }else if (boxplot.aggr %in% c("c","unlist")){
                  original = eft[[var.vname]]$tspftdbh[,,y.sel,,dbh.bp,n.pft]
                  concat   = apply( X      = expand.grid( dimnames(original)[[2]]
                                                        , dimnames(original)[[3]]
                                                        )#end expand.grid
                                  , MARGIN = 1
                                  , FUN    = paste
                                  , collapse = "+"
                                  )#end apply
                  every    = apply( X      = original
                                  , MARGIN = c(1,4,5)
                                  , FUN    = unlist
                                  )#end apply
                  every    = aperm(a=every,perm=c(2,1,3,4))
                  dimnames(every)[[2]] = concat
               }else{
                  every = apply( X      = eft[[var.vname]]$tspftdbh[,,y.sel,,dbh.bp,n.pft]
                               , MARGIN = c(1,2,4,5)
                               , FUN    = get(boxplot.aggr)
                               , na.rm  = TRUE
                               )#end apply
               }#end if
               ylimit.season = pretty.xylim( u       = every[l.sel,,-n.season,]
                                           , fracexp = 0.0
                                           , is.log  = var.plog.dbh
                                           )#end pretty.xylim
               ylimit.year   = pretty.xylim( u       = every[l.sel,, n.season,]
                                           , fracexp = 0.0
                                           , is.log  = var.plog.dbh
                                           )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #----- Load variable and set y axis. ---------------------------------------#
               this   = every[sel,,-n.season,]
               ylimit = ylimit.season
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Seasonal means) ","\n",out.desc,sep="")
               lex     = lex     = desc.unit(desc="DBH class",unit=untab$cm)
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the box plot colour.                                           #
               #---------------------------------------------------------------------------#
               bp.colour = rep(c("transparent",now$colour,"transparent"),times=n.dbh.bp)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$boxdbh.season
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-boxdbh_season.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wsize$width,height=wsize$height
                               ,pointsize=ptsz,paper=wsize$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=wsize$width,height=wsize$height
                        ,pointsize=ptsz,paper=wsize$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(par.user)
                  par(oma = c(0.2,3,4.5,0))
                  layout(mat    = rbind(lo.season$mat.off,rep(1,times=lo.season$ncol))
                        ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , fill    = now$colour
                         , border  = foreground
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (e in 1:n.season.mp){

                     #---------------------------------------------------------------------#
                     #     Find out where is this box going, and set up axes and margins.  #
                     #---------------------------------------------------------------------#
                     left    = (e %% lo.season$ncol) == 1
                     right   = (e %% lo.season$ncol) == 0
                     top     = e <= lo.season$ncol
                     bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(season.desc[e],sep="")
                     #---------------------------------------------------------------------#



                     #----- Create the temporary box plot list. ---------------------------#
                     ss.this = this[,,e,]
                     empty   = array( data     = NA
                                    , dim      = dim(ss.this)[-1]
                                    , dimnames = dimnames(ss.this)[-1]
                                    )#end array
                     ss.this = abind(aaa=empty,ss.this,zzz=empty,along=1)
                     ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                     bp.plot = split(x=ss.this,f=list(ai[,1],ai[,3]))
                     #---------------------------------------------------------------------#



                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n"
                                ,log=var.xylog.dbh)
                     if (bottom) axis(side=1,at=xat,labels=dbh.bp.key)
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")

                     #----- Add the box plot, without the x axis. -------------------------#
                     boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                            ,yaxt="n")
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.season.mp)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 1/6
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff.e)
                  # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the box plot by DBH (annual means).                             #
               #---------------------------------------------------------------------------#
               xlimit = c(0,n.dbh.bp*n.lax) + 0.5
               xat    = seq(from=1 + 0.5*(n.lax-1), to=n.dbh.bp*n.lax, by=n.lax)
               xgrid  = seq(from=0.5, to=(n.dbh.bp+1)*n.lax - 0.5*(n.lax-1), by=n.lax)
               #---------------------------------------------------------------------------#



               #----- Load variable and set y axis. ---------------------------------------#
               this   = every[sel,, n.season,]
               ylimit = ylimit.year
               #---------------------------------------------------------------------------#



               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means)","\n",out.desc,sep="")
               lex     = desc.unit(desc="DBH class",unit=untab$cm)
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Set up the box plot colour.                                           #
               #---------------------------------------------------------------------------#
               bp.colour = rep(c("transparent",now$colour,"transparent"),times=n.dbh.bp)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$boxdbh.year
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix,"-boxdbh_year."
                                 ,outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Create the temporary box plot list. ------------------------------#
                  empty   = array( data     = NA
                                 , dim      = dim(this)[-1]
                                 , dimnames = dimnames(this)[-1]
                                 )#end array
                  ss.this = abind(aaa=empty,this,zzz=empty,along=1)
                  ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                  bp.plot = split(x=ss.this,f=list(ai[,1],ai[,3]))
                  #------------------------------------------------------------------------#


                  #----- Split into two plots. --------------------------------------------#
                  par(par.user)
                  layout(mat=matrix(c(2,1),nrow=2),heights=c(5,1))
                  #------------------------------------------------------------------------#



                  #----- Add the legend. --------------------------------------------------#
                  par(mar=c(0.2,4.4,0.1,2.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = inset
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , fill    = now$colour
                         , border  = foreground
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=c(5.1,4.4,4.1,2.1))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n",log=var.xylog.dbh)
                  axis(side=1,at=xat,labels=dbh.bp.key)
                  axis(side=2)
                  box()
                  title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.ptsz)
                  if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")

                  #----- Add the box plot, without the x axis. ----------------------------#
                  boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                         ,yaxt="n")
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (var.plt & is.dbh)
            #..............................................................................#
            #..............................................................................#






            #..............................................................................#
            #..............................................................................#
            #     Boxplots by PFT, DBH, and season...                                      #
            #------------------------------------------------------------------------------#
            if (var.plt & is.dbh){
               #---------------------------------------------------------------------------#
               #     Load all simulations for y axis.                                      #
               #---------------------------------------------------------------------------#
               if (n.real == 1){
                  every = eft[[var.vname]]$tspftdbh[,1,y.sel,,dbh.bp,pft.mp]
               }else if (boxplot.aggr %in% c("c","unlist")){
                  original = eft[[var.vname]]$tspftdbh[,,y.sel,,dbh.bp,pft.mp]
                  concat   = apply( X      = expand.grid( dimnames(original)[[2]]
                                                        , dimnames(original)[[3]]
                                                        )#end expand.grid
                                  , MARGIN = 1
                                  , FUN    = paste
                                  , collapse = "+"
                                  )#end apply
                  every    = apply( X      = original
                                  , MARGIN = c(1,4,5,6)
                                  , FUN    = unlist
                                  )#end apply
                  every    = aperm(a=every,perm=c(2,1,3,4,5))
                  dimnames(every)[[2]] = concat
               }else{
                  every = apply( X      = eft[[var.vname]]$tspftdbh[,,y.sel,,dbh.bp,pft.mp]
                               , MARGIN = c(1,2,4,5,6)
                               , FUN    = get(boxplot.aggr)
                               , na.rm  = TRUE
                               )#end apply
               }#end if
               ylimit.season = pretty.xylim( u       = every[l.sel,,-n.season,,]
                                           , fracexp = 0.0
                                           , is.log  = var.plog.dbh
                                           )#end pretty.xylim
               ylimit.year   = pretty.xylim( u       = every[l.sel,, n.season,,]
                                           , fracexp = 0.0
                                           , is.log  = var.plog.dbh
                                           )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over each season.                                                #
               #---------------------------------------------------------------------------#
               for (e in 1:n.season){

                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                  #      Plot the box plot by PFT and by season.                           #
                  #------------------------------------------------------------------------#
                  xlimit = c(0,n.dbh.bp*n.lax) + 0.5
                  xat    = seq(from=1 + 0.5*(n.lax-1), to=n.dbh.bp*n.lax, by=n.lax)
                  xgrid  = seq(from=0.5,to=(n.dbh.bp+1)*n.lax-0.5*(n.lax-1),by=n.lax)
                  #------------------------------------------------------------------------#




                  #----- Load variable and set y axis. ------------------------------------#
                  this   = every[sel,,e,,]
                  if (e == n.season){
                     ylimit = ylimit.year
                  }else{
                     ylimit = ylimit.season
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Set the title. ---------------------------------------------------#
                  letitre = paste(var.desc," - ","Season: ",season.desc[e]
                                 ,"\n",out.desc,sep="")
                  lex     = desc.unit(desc="DBH class",unit=untab$cm)
                  ley     = desc.unit(desc=var.desc,unit=var.unit)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Set up the box plot colour.                                        #
                  #------------------------------------------------------------------------#
                  bp.colour = rep(c("transparent",now$colour,"transparent"),times=n.dbh.bp)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all formats.                                            #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Get the path. -------------------------------------------------#
                     now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$boxpftdbh[e]
                     #---------------------------------------------------------------------#


                     #----- Open file or display. -----------------------------------------#
                     fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                    ,"-boxpftdbh-year.",outform[o],sep="")
                     if (outform[o] == "x11"){
                        X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=wsize$width*depth
                           ,height=wsize$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=wsize$width,height=wsize$height
                                  ,pointsize=ptsz,paper=wsize$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE,width=wsize$width
                           ,height=wsize$height,pointsize=ptsz,paper=wsize$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Split the window into several smaller windows.  Add a bottom    #
                     # row to fit the legend.                                              #
                     #---------------------------------------------------------------------#
                     par.orig = par(no.readonly = TRUE)
                     mar.orig = par.orig$mar
                     par(par.user)
                     par(oma = c(0.2,3,4.5,0))
                     layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                           ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                           )#end layout
                     #---------------------------------------------------------------------#



                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.2,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = if (any(now$parse)){
                                           parse(text=now$legend)
                                        }else{
                                           now$legend
                                        }#end if
                            , fill    = now$colour
                            , border  = foreground
                            , bg      = background
                            , ncol    = min(3,pretty.box(n.now)$ncol)
                            , title   = expression(bold("Simulation"))
                            , cex     = cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Loop over all seasons, and plot the bar plots.                  #
                     #---------------------------------------------------------------------#
                     for (f in 1:n.pft.mp){
                        #------------------------------------------------------------------#
                        #       Find out where is this box going, and set up axes and      #
                        #  margins.                                                        #
                        #------------------------------------------------------------------#
                        left    = (f %% lo.pft$ncol) == 1
                        right   = (f %% lo.pft$ncol) == 0
                        top     = f <= lo.pft$ncol
                        bottom  = f > (lo.pft$nrow - 1) * lo.pft$ncol
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                        #------------------------------------------------------------------#


                        #----- Set up the title for each plot. ----------------------------#
                        lesub = paste(pft.desc[pft.mp[f]],sep="")
                        #------------------------------------------------------------------#



                        #----- Create the temporary box plot list. ------------------------#
                        ss.this = this[,,,f]
                        empty   = array( data     = NA
                                       , dim      = dim(ss.this)[-1]
                                       , dimnames = dimnames(ss.this)[-1]
                                       )#end array
                        ss.this = abind(aaa=empty,ss.this,zzz=empty,along=1)
                        ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                        bp.plot = split(x=ss.this,f=list(ai[,1],ai[,3]))
                        #------------------------------------------------------------------#



                        #----- Plot window and grid. --------------------------------------#
                        par(mar=mar.now)
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n"
                                   ,log=var.xylog.dbh)
                        if (bottom) axis(side=1,at=xat,labels=dbh.bp.key)
                        if (left  ) axis(side=2)
                        box()
                        title(main=lesub,cex.main=cex.ptsz,line=0.5)
                        if (plotgrid){
                           abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                        }#end if

                        #----- Add the box plot, without the x axis. ----------------------#
                        boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE
                               ,show.names=FALSE,yaxt="n")
                        #------------------------------------------------------------------#
                     }#end for (f in 1:n.pft.mp)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the global title.                                          #
                     #---------------------------------------------------------------------#
                     gtitle( main      = letitre
                           , xlab      = lex
                           , ylab      = ley
                           , off.xlab  = 1/6
                           , line.xlab = 4.1
                           , line.ylab = 2.6
                           , cex.main  = 1.1*cex.ptsz
                           )#end gtitle
                     # par(las=0)
                     # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff.e)
                     # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                     # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     clean.tmp()
                     #---------------------------------------------------------------------#
                  }#end for (o in 1:nout)
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               }#end for (e in 1:n.season)
               #---------------------------------------------------------------------------#
            }#end if (var.plt && is.dbh)
            #------------------------------------------------------------------------------#
         }#end for
         #=================================================================================#
         #=================================================================================#
         #.................................................................................#
         #.................................................................................#






         #.................................................................................#
         #.................................................................................#
         #      Plot the bar plots.                                                        #
         #---------------------------------------------------------------------------------#
         cat  ("         ~ Plotting the bar plots...","\n")
         for (v in loop.barplot){
            #----- Copy variable info. ----------------------------------------------------#
            var.vname  = scen.barplot$vname [v]
            var.desc   = scen.barplot$desc  [v]
            var.unit   = scen.barplot$unit  [v]
            is.pft     = scen.barplot$pftvar[v]
            is.dbh     = scen.barplot$dbhvar[v]
            var.plog   = scen.barplot$plog  [v]
            var.plt    = scen.barplot$plt   [v]
            if (var.plog){
               var.xylog  = "y"
            }else{
               var.xylog  = ""
            }#end if
            #------------------------------------------------------------------------------#


            if (var.plt){
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the bar plot by PFT, DBH class, and season.                     #
               #---------------------------------------------------------------------------#
               every = eft[[var.vname]]$tspftdbh[,,y.sel,,,pft.bp,drop=FALSE]
               every = apply( X      = every
                            , MARGIN = c(1,4,5,6)
                            , FUN    = get(barplot.aggr)
                            , na.rm  = TRUE
                            )#end apply
               every = apply( X = every, MARGIN = c(1,2,3)  , FUN = cumsum)
               every = aperm( a = every, perm   = c(2,4,1,3)              )
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Find the y limits:                                                   #
               #---------------------------------------------------------------------------#
               yrange.season = pretty.xylim( u       = c(0,every[l.sel,,,-n.season])
                                           , fracexp =  0.0
                                           , is.log  = FALSE
                                           )#end pretty.xylim
               ylimit.season = pretty.xylim( u       = c(0,every[l.sel,,,-n.season])
                                           , fracexp = -0.06
                                           , is.log  = FALSE
                                           )#end pretty.xylim
               yrange.year   = pretty.xylim( u       = c(0,every[l.sel,,, n.season])
                                           , fracexp =  0.0
                                           , is.log  = FALSE
                                           )#end pretty.xylim
               ylimit.year   = pretty.xylim( u       = c(0,every[l.sel,,, n.season])
                                           , fracexp = -0.06
                                           , is.log  = FALSE
                                           )#end pretty.xylim
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Stack the PFT data.                                                  #
               #---------------------------------------------------------------------------#
               this                = every[sel,,,]
               dimnames(this)[[1]] = now$key
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Combine simulation and size into one dimension.                      #
               #---------------------------------------------------------------------------#
               dth  = dim(this)
               dnth = list( paste(rep(dimnames(this)[[1]],times=dth[2])
                                 ,rep(dimnames(this)[[2]],each =dth[1])
                                 ,sep=".")
                          , dimnames(this)[[3]]
                          , dimnames(this)[[4]]
                          )#end list
               this = array(data= this,dim=c(dth[1]*dth[2],dth[3],dth[4]),dimnames=dnth)
               #---------------------------------------------------------------------------#


               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Seasonal means)","\n",out.desc,sep="")
               lex     = desc.unit(desc="DBH class",unit=untab$cm)
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    Find the x dimensions.                                                 #
               #---------------------------------------------------------------------------#
               xleft  = ( rep(sequence(n.now),times=n.dbh)
                        + rep((n.now+1)*(sequence(n.dbh)-1),each=n.now) ) - 0.9
               xright = xleft + 0.8
               xat    = seq(from=0.5*n.now,by=n.now+1,length.out=n.dbh)
               xgrid  = seq(from=0,by=n.now+1,length.out=n.dbh+1)-0.5
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find out the colour filling.                                          #
               #---------------------------------------------------------------------------#
               fill.col   = array( data     = rep(pft.colour[pft.bp],each=dim(this)[1])
                                 , dim      = dim(this)[1:2]
                                 , dimnames = dimnames(this)[1:2]
                                 )#end array
               refer.col  = array( data     = rep(now$colour,times=dim(this)[2])
                                 , dim      = dim(this)[1:2]
                                 , dimnames = dimnames(this)[1:2]
                                 )#end array
               xleft      = array( data     = xleft
                                 , dim      = dim(this)[1:2]
                                 , dimnames = dimnames(this)[1:2]
                                 )#end array
               xright     = array( data     = xright
                                 , dim      = dim(this)[1:2]
                                 , dimnames = dimnames(this)[1:2]
                                 )#end array
               #------- Fix the y for stacked. --------------------------------------------#
               ybottom       = this
               ybottom[,-1,] = ybottom[,-dim(ybottom)[2],]
               ybottom[ ,1,] = 0
               ytop          = this
               #---------------------------------------------------------------------------#




               #----- Find the limits for the plot. ---------------------------------------#
               xlimit = range(c(xat,xgrid))
               yrange = yrange.season
               ylimit = ylimit.season
               #---------------------------------------------------------------------------#



               #------- Fix the y for stacked. --------------------------------------------#
               lbottom = yrange[1] - 0.06 * diff(yrange) + 0 * ybottom
               ltop    = yrange[1] - 0.02 * diff(yrange) + 0 * ybottom
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$barplot.season
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-barplot_season.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wsize$width,height=wsize$height
                               ,pointsize=ptsz,paper=wsize$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=wsize$width,height=wsize$height
                        ,pointsize=ptsz,paper=wsize$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = c(5.1,4.4,4.1,2.1)
                  par(par.user)
                  par(oma = c(0.2,3,4.5,0))

                  emat = matrix( data  = c( rep(t(lo.season$mat.off2),each=2)
                                          , rep(1,times=lo.season$ncol)
                                          , rep(2,times=lo.season$ncol)
                                          )#end c
                               , ncol  = 2*lo.season$ncol
                               , nrow  = 1+lo.season$nrow
                               , byrow = TRUE
                               )#end matrix

                  layout(mat     = emat
                        ,heights = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot simulation legend. ------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , fill    = now$colour
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #----- Plot PFT legend. -------------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = pft.desc  [pft.bp]
                         , fill    = pft.colour[pft.bp]
                         , border  = foreground
                         , bg      = background
                         , ncol    = min(2,pretty.box(n.pft.bp)$ncol)
                         , title   = expression(bold("Plant Functional Type"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Loop over all seasons, and plot the bar plots.                     #
                  #------------------------------------------------------------------------#
                  for (e in 1:n.season.mp){

                     #---------------------------------------------------------------------#
                     #    Find out where is this box going, and set up axes and margins.   #
                     #---------------------------------------------------------------------#
                     left    = (e %% lo.season$ncol) == 1
                     right   = (e %% lo.season$ncol) == 0
                     top     = e <= lo.season$ncol
                     bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     #---------------------------------------------------------------------#


                     #----- Set up the title for each plot. -------------------------------#
                     lesub = paste(season.desc[e],sep="")
                     #---------------------------------------------------------------------#



                     #----- Plot window and grid. -----------------------------------------#
                     par(mar=mar.now)
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                     if (bottom) axis(side=1,at=xat,labels=dbh.key)
                     if (left  ) axis(side=2)
                     box()
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                     #---------------------------------------------------------------------#



                     #----- Add the bar plot using rect, so we can change the line width. -#
                     rect( xleft   = xleft
                         , ybottom = ybottom [,,e]
                         , xright  = xright
                         , ytop    = ytop    [,,e]
                         , density = -1
                         , col     = fill.col
                         , border  = foreground
                         , lwd     = barplot.lwd
                         )#end rect
                     #---------------------------------------------------------------------#



                     #----- Add the simulation colour beneath the bar plot. ---------------#
                     rect( xleft   = xleft
                         , ybottom = lbottom [,,e]
                         , xright  = xright
                         , ytop    = ltop    [,,e]
                         , density = -1
                         , col     = refer.col
                         , border  = refer.col
                         , lwd     = barplot.lwd
                         )#end rect
                     #---------------------------------------------------------------------#
                  }#end for (e in 1:n.season.mp)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 1/6
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff.e)
                  # mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # mtext(text=letitre,side=3,outer=TRUE,cex=1.1*cex.ptsz,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               #      Plot the bar plot by PFT, DBH class, using annual means.             #
               #---------------------------------------------------------------------------#


               #----- Set the title. ------------------------------------------------------#
               letitre = paste(var.desc," (Annual means)","\n",out.desc,sep="")
               lex     = desc.unit(desc="DBH class",unit=untab$cm)
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#




               #----- Find the limits for the plot. ---------------------------------------#
               xlimit = range(c(xat,xgrid))
               yrange = yrange.year
               ylimit = ylimit.year
               #---------------------------------------------------------------------------#



               #------- Fix the y for stacked. --------------------------------------------#
               lbottom = yrange[1] - 0.06 * diff(yrange) + 0 * ybottom
               ltop    = yrange[1] - 0.02 * diff(yrange) + 0 * ybottom
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path. ----------------------------------------------------#
                  now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$barplot.year
                  #------------------------------------------------------------------------#


                  #----- Open file or display. --------------------------------------------#
                  fichier = paste(now.outpath,"/",var.vname,"-",out.suffix
                                 ,"-barplot_year.",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.  Add a bottom row   #
                  # to fit the legend.                                                     #
                  #------------------------------------------------------------------------#
                  par.orig = par(no.readonly = TRUE)
                  emat = matrix(c(3,3,1,2),ncol=2,nrow=2,byrow=TRUE)
                  par(par.user)
                  layout(mat=emat,heights=c(5,1),widths=c(3,2))
                  #------------------------------------------------------------------------#



                  #----- Plot simulation legend. ------------------------------------------#
                  par(mar=c(0.2,4.4,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = if (any(now$parse)){
                                        parse(text=now$legend)
                                     }else{
                                        now$legend
                                     }#end if
                         , fill    = now$colour
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.now)$ncol)
                         , title   = expression(bold("Simulation"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #----- Plot PFT legend. -------------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,2.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = pft.desc  [pft.bp]
                         , fill    = pft.colour[pft.bp]
                         , border  = foreground
                         , bg      = background
                         , ncol    = min(3,pretty.box(n.pft.bp)$ncol)
                         , title   = expression(bold("Plant Functional Type"))
                         , cex     = cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#




                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=mar.orig)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  axis(side=1,at=xat,labels=dbh.key)
                  axis(side=2)
                  box()
                  title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.ptsz)
                  if (plotgrid) abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                  #------------------------------------------------------------------------#



                  #----- Add the bar plot using rect, so we can change the line width. ----#
                  rect( xleft   = xleft
                      , ybottom = ybottom[,,n.season]
                      , xright  = xright
                      , ytop    = ytop   [,,n.season]
                      , density = -1
                      , col     = fill.col
                      , border  = foreground
                      , lwd     = barplot.lwd
                      )#end rect
                  #------------------------------------------------------------------------#



                  #----- Add the simulation colour beneath the bars. ----------------------#
                  rect( xleft   = xleft
                      , ybottom = lbottom[,,n.season]
                      , xright  = xright
                      , ytop    = ltop   [,,n.season]
                      , density = -1
                      , col     = refer.col
                      , border  = refer.col
                      , lwd     = barplot.lwd
                      )#end rect
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (var.plt)
            #------------------------------------------------------------------------------#
         }#end for (v in 1:nscen.barplot)
         #.................................................................................#
         #.................................................................................#






         #.................................................................................#
         #.................................................................................#
         #      Plot the parameter space.                                                  #
         #---------------------------------------------------------------------------------#
         cat  ("         ~ Plotting the parameter space...","\n")
         #---------------------------------------------------------------------------------#
         #     Loop over Y variables.                                                      #
         #---------------------------------------------------------------------------------#
         for (y in loop.y){
            y.vname    = scen.xyz$yvar$vname   [y]
            y.desc     = scen.xyz$yvar$desc    [y]
            y.unit     = scen.xyz$yvar$unit    [y]
            y.plog     = scen.xyz$yvar$plog    [y]
            y.plog.dbh = scen.xyz$yvar$plog.dbh[y]
            y.leg      = scen.xyz$yvar$leg     [y]
            y.pft      = scen.xyz$yvar$pftvar  [y]
            y.dbh      = scen.xyz$yvar$dbhvar  [y]

            cat  ("           = Y variable: ",y.desc,"...","\n")

            #----- Time series. -----------------------------------------------------------#
            y.ts      = eft[[y.vname]]$ts      [sel,,y.sel,season.mp]
            #------------------------------------------------------------------------------#


            #----- Make sure the PFT version is PFT-dependent. ----------------------------#
            if (y.pft){
               y.tspft = eft[[y.vname]]$tspft[sel,,y.sel,n.season,pft.mp]
            }else{
               y.tspft    = eft[[y.vname]]$ts      [sel,,y.sel,n.season]
               y.dim      = c(dim(y.tspft),n.pft.mp)
               y.dimnames = c(dimnames(y.tspft),list(pft.mp.key))
               y.tspft    = array(NA,dim=y.dim,dimnames=y.dimnames)
               for (f in 1:n.pft.mp){
                  y.tspft[,,,f] = eft[[y.vname]]$ts[sel,,y.sel,n.season]
               }#end for
            }#end if
      
      
      
            #------------------------------------------------------------------------------#


            #----- Make sure that the DBH version is DBH-dependent. -----------------------#
            if (y.dbh){
               y.tsdbh = eft[[y.vname]]$tspftdbh[sel,,y.sel,n.season,dbh.mp,n.pft]
            }else{
               y.tsdbh    = eft[[y.vname]]$ts      [sel,,y.sel,n.season]
               y.dim      = c(dim(y.tsdbh),n.dbh.mp)
               y.dimnames = c(dimnames(y.tsdbh),list(dbh.mp.key))
               y.tsdbh    = array(NA,dim=y.dim,dimnames=y.dimnames)
               for (d in 1:n.dbh.mp){
                  y.tsdbh[,,,d] = eft[[y.vname]]$ts[sel,,y.sel,n.season]
               }#end for
            }#end if
            #------------------------------------------------------------------------------#


            #----- Make sure that the DBH+PFT version has both DBH and PFT dimensions. ----#
            if (y.dbh){
               y.tsall = eft[[y.vname]]$tspftdbh[sel,,y.sel,,,pft.mp]
            }else{
               y.tsall    = eft[[y.vname]]$ts      [sel,,y.sel,]
               y.dim      = c(dim(y.tsall),n.dbh,n.pft.mp)
               y.dimnames = c(dimnames(y.tsall),list(dbh.key),list(pft.mp.key))
               y.tsall = array(NA,dim=y.dim,dimnames=y.dimnames)
               for (d in 1:n.dbh){
                  for (f in 1:n.pft.mp){
                     y.tsall[,,,,d,f] = eft[[y.vname]]$ts[sel,,y.sel,]
                  }#end for
               }#end for
            }#end if
            #------------------------------------------------------------------------------#


            #----- Create indices to split the arrays into lists. -------------------------#
            ai.ts    = arrayInd(sequence(length(y.ts   )),.dim=dim(y.ts   ))
            ai.tspft = arrayInd(sequence(length(y.tspft)),.dim=dim(y.tspft))
            ai.tsdbh = arrayInd(sequence(length(y.tsdbh)),.dim=dim(y.tsdbh))
            #------------------------------------------------------------------------------#


            #----- Find the format of the points for the plots. ---------------------------#
            pch.ts    = array( data = now$pch[arrayInd(1:length(y.ts   )
                                                      ,.dim=dim(y.ts   ))[,1]]
                             , dim  = dim(y.ts)
                             )#end array
            pch.tspft = array( data = now$pch[arrayInd( 1:length(y.tspft)
                                                      , .dim=dim(y.tspft))[,1]]
                             , dim  = dim(y.tspft)
                             )#end array
            pch.tsdbh = array( data = now$pch[arrayInd( 1:length(y.tsdbh)
                                                      , .dim=dim(y.tsdbh))[,1]]
                             , dim  = dim(y.tsdbh)
                             )#end array
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Loop over X variables.                                                   #
            #------------------------------------------------------------------------------#
            for (x in loop.x){
               x.vname    = scen.xyz$xvar$vname   [x]
               x.desc     = scen.xyz$xvar$desc    [x]
               x.unit     = scen.xyz$xvar$unit    [x]
               x.plog     = scen.xyz$xvar$plog    [x]
               x.plog.dbh = scen.xyz$xvar$plog.dbh[x]
               x.leg      = scen.xyz$xvar$leg     [x]
               x.pft      = scen.xyz$xvar$pftvar  [x]
               x.dbh      = scen.xyz$xvar$dbhvar  [x]


               #----- Time series. --------------------------------------------------------#
               x.ts       = eft[[x.vname]]$ts     [sel,,y.sel,season.mp]
               #---------------------------------------------------------------------------#


               #----- Make sure the PFT version is PFT-dependent. -------------------------#
               if (x.pft){
                  x.tspft = eft[[x.vname]]$tspft[sel,,y.sel,n.season,pft.mp]
               }else{
                  x.tspft = array(NA,dim(y.tspft),dimnames=dimnames(y.tspft))
                  for (f in 1:n.pft.mp){
                     x.tspft[,,,f] = eft[[x.vname]]$ts[sel,,y.sel,n.season]
                  }#end for
               }#end if
               #---------------------------------------------------------------------------#


               #----- Make sure the DBH version is DBH-dependent. -------------------------#
               if (x.dbh){
                  x.tsdbh = eft[[x.vname]]$tspftdbh[sel,,y.sel,n.season,dbh.mp,n.pft]
               }else{
                  x.tsdbh = array(NA,dim(y.tsdbh),dimnames=dimnames(y.tsdbh))
                  for (d in 1:n.dbh.mp){
                     x.tsdbh[,,,d] = eft[[x.vname]]$ts[sel,,y.sel,n.season]
                  }#end for
               }#end if
               #---------------------------------------------------------------------------#


               #----- Make sure the DBH+PFT version depends on both DBH and PFT. ----------#
               if (x.dbh){
                  x.tsall = eft[[x.vname]]$tspftdbh[sel,,y.sel,,,pft.mp]
               }else{
                  x.tsall = array(NA,dim(y.tsall),dimnames=dimnames(y.tsall))
                  for (d in 1:n.dbh){
                     for (f in 1:n.pft.mp){
                        x.tsall[,,,,d,f] = eft[[x.vname]]$ts[sel,,y.sel,]
                     }#end for
                  }#end for
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Plot only when x and y are different.                                 #
               #---------------------------------------------------------------------------#
               if (x.vname != y.vname){
                  #----- Build the log info. ----------------------------------------------#
                  xy.plog=""
                  if (x.plog) xy.plog=paste(xy.plog,"x",sep="")
                  if (y.plog) xy.plog=paste(xy.plog,"y",sep="")
                  xy.plog.dbh=""
                  if (x.plog.dbh) xy.plog.dbh=paste(xy.plog.dbh,"x",sep="")
                  if (y.plog.dbh) xy.plog.dbh=paste(xy.plog.dbh,"y",sep="")
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over Z variables.                                             #
                  #------------------------------------------------------------------------#
                  for (z in loop.z){
                     z.vname    = scen.xyz$zvar$vname         [z]
                     z.desc     = scen.xyz$zvar$desc          [z]
                     z.unit     = scen.xyz$zvar$unit          [z]
                     z.plog     = scen.xyz$zvar$plog          [z]
                     z.plog.dbh = scen.xyz$zvar$plog.dbh      [z]
                     z.pft      = scen.xyz$zvar$pftvar        [z]
                     z.dbh      = scen.xyz$zvar$dbhvar        [z]
                     z.cscheme  = get(scen.xyz$zvar$col.scheme[z])


                     #---------------------------------------------------------------------#
                     #     Plot only when x, y, and z are different.                       #
                     #---------------------------------------------------------------------#
                     if (! z.vname %in% c(x.vname,y.vname)){
                        #----- Time series. -----------------------------------------------#
                        z.ts       = eft[[z.vname]]$ts   [sel,,y.sel,season.mp]
                        #------------------------------------------------------------------#


                        #----- Make sure the PFT version is PFT-dependent. ----------------#
                        if (z.pft){
                           z.tspft = eft[[z.vname]]$tspft[sel,,y.sel,n.season,pft.mp]
                        }else{
                           z.tspft = array(NA,dim=dim(y.tspft),dimnames=dimnames(y.tspft))
                           for (f in 1:n.pft.mp){
                              z.tspft[,,,f] = eft[[z.vname]]$ts[sel,,y.sel,n.season]
                           }#end z.tspft
                        }#end if (z.pft)
                        #------------------------------------------------------------------#


                        #----- Make sure the DBH version is DBH-dependent. ----------------#
                        if (z.dbh){
                           z.tsdbh = eft[[z.vname]]$tspftdbh
                           z.tsdbh = z.tsdbh[sel,,y.sel,n.season,dbh.mp,n.pft]
                        }else{
                           z.tsdbh = array(NA,dim=dim(y.tsdbh),dimnames=dimnames(y.tsdbh))
                           for (d in 1:n.dbh.mp){
                              z.tsdbh[,,,d] = eft[[z.vname]]$ts[sel,,y.sel,n.season]
                           }#end z.tspft
                        }#end if (z.dbh)
                        #------------------------------------------------------------------#


                        #----- Make sure the DBH+PFT version depends on both DBH and PFT. -#
                        if (z.dbh){
                           z.tsall = eft[[z.vname]]$tspftdbh[sel,,y.sel,,,pft.mp]
                        }else{
                           z.tsall = array(NA,dim(y.tsall),dimnames=dimnames(y.tsall))
                           for (d in 1:n.dbh){
                              for (f in 1:n.pft.mp){
                                 z.tsall[,,,,d,f] = eft[[z.vname]]$ts[sel,,y.sel,]
                              }#end for
                           }#end for
                        }#end if
                        #------------------------------------------------------------------#



                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #     Split the arrays into lists.                                 #
                        #------------------------------------------------------------------#
                        x.list    = split(x=x.ts  ,f=ai.ts[,4])
                        y.list    = split(x=y.ts  ,f=ai.ts[,4])
                        z.list    = split(x=z.ts  ,f=ai.ts[,4])
                        pch.list  = split(x=pch.ts,f=ai.ts[,4])
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Fix the titles.                                             #
                        #------------------------------------------------------------------#
                        letitre   = paste("Seasonal means - ",z.desc,"\n",out.desc,sep="")
                        lex       = desc.unit(desc=x.desc,unit=x.unit)
                        ley       = desc.unit(desc=y.desc,unit=y.unit)
                        lacle     = desc.unit(desc=NULL  ,unit=z.unit)
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Plot the parameter space by season.                          #
                        #------------------------------------------------------------------#
                        for (o in 1:nout){
                           #----- Get the path. -------------------------------------------#
                           now.scenario = outpath[[o]]$global[[g]]$scenario[[s]]
                           now.outpath  = now.scenario$xyz.season[y]
                           #---------------------------------------------------------------#



                           #----- Open file or display. -----------------------------------#
                           fichier = paste(now.outpath,"/y-",y.vname,"_x-",x.vname
                                          ,"_z-",z.vname,"-",out.suffix,"-season."
                                          ,outform[o],sep="")
                           if (outform[o] == "x11"){
                              X11(width=size$width,height=size$height,pointsize=ptsz)
                           }else if(outform[o] == "png"){
                              png(filename=fichier,width=size$width*depth
                                 ,height=size$height*depth,pointsize=ptsz,res=depth)
                           }else if(outform[o] == "eps"){
                              postscript(file=fichier,width=size$width,height=size$height
                                        ,pointsize=ptsz,paper=size$paper)
                           }else if(outform[o] == "pdf"){
                              pdf(file=fichier,onefile=FALSE,width=size$width
                                 ,height=size$height,pointsize=ptsz,paper=size$paper)
                           }#end if
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #      Plot the boxes.                                          #
                           #---------------------------------------------------------------#
                           leg.ncol    = min(3,pretty.box(n.now)$ncol)
                           leg.letitre = expression(bold("Simulation"))
                           leg.legend  = if (any(now$parse)){
                                            parse(text=now$legend)
                                         }else{
                                            now$legend
                                         }#end if
                           par(par.user)
                           xyz.plot( x              = x.list
                                   , y              = y.list
                                   , z              = z.list
                                   , fixed.xlim     = TRUE
                                   , fixed.ylim     = TRUE
                                   , xy.log         = xy.plog
                                   , pch            = pch.list
                                   , cex            = 1.2
                                   , nlevels        = n.colourbar
                                   , colour.palette = z.cscheme
                                   , xyz.sub        = season.mp.desc
                                   , xyz.title      = list( main      = letitre
                                                          , xlab      = lex
                                                          , ylab      = ley
                                                          , line.xlab = 4.1
                                                          , line.ylab = 2.6
                                                          , cex.main  = 1.1 * cex.ptsz
                                                          )#end list
                                   , xyz.more       = list(grid=list(col=grid.colour
                                                                    ,lty="solid"))
                                   , key.log        = z.plog
                                   , key.title      = list(main=lacle,cex.main=.8*cex.main)
                                   , xyz.legend     = list( x      = "bottom"
                                                          , inset  = 0.0
                                                          , legend = leg.legend
                                                          , pch    = now$pch
                                                          , col    = grey.fg
                                                          , border = foreground
                                                          , bg     = background
                                                          , ncol   = leg.ncol
                                                          , title  = leg.letitre
                                                          , cex    = 1.0
                                                          , xpd    = TRUE
                                                          )#end legend
                                   )#end xyz.plot
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           clean.tmp()
                           #---------------------------------------------------------------#
                        }#end for
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#











                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #     XYZ plots by PFT (annual means).                             #
                        #------------------------------------------------------------------#



                        #----- Split the arrays into lists. -------------------------------#
                        x.list    = split(x=x.tspft  ,f=ai.tspft[,4])
                        y.list    = split(x=y.tspft  ,f=ai.tspft[,4])
                        z.list    = split(x=z.tspft  ,f=ai.tspft[,4])
                        pch.list  = split(x=pch.tspft,f=ai.tspft[,4])
                        #------------------------------------------------------------------#



                        #----- Fix the titles. --------------------------------------------#
                        letitre   = paste("Annual means - ",z.desc,"\n",out.desc,sep="")
                        lex       = desc.unit(desc=x.desc,unit=x.unit)
                        ley       = desc.unit(desc=y.desc,unit=y.unit)
                        lacle     = desc.unit(desc=NULL  ,unit=z.unit)
                        #------------------------------------------------------------------#




                        #----- Plot the parameter space by PFT. ---------------------------#
                        for (o in 1:nout){
                           #----- Get the path. -------------------------------------------#
                           now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$xyz.pft[y]
                           #---------------------------------------------------------------#

                           #----- Open file or display. -----------------------------------#
                           fichier = paste(now.outpath,"/y-",y.vname,"_x-",x.vname
                                          ,"_z-",z.vname,"-",out.suffix,"-pft."
                                          ,outform[o],sep="")
                           if (outform[o] == "x11"){
                              X11(width=size$width,height=size$height,pointsize=ptsz)
                           }else if(outform[o] == "png"){
                              png(filename=fichier,width=size$width*depth
                                 ,height=size$height*depth,pointsize=ptsz,res=depth)
                           }else if(outform[o] == "eps"){
                              postscript(file=fichier,width=size$width,height=size$height
                                        ,pointsize=ptsz,paper=size$paper)
                           }else if(outform[o] == "pdf"){
                              pdf(file=fichier,onefile=FALSE,width=size$width
                                 ,height=size$height,pointsize=ptsz,paper=size$paper)
                           }#end if
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #      Plot the boxes.                                          #
                           #---------------------------------------------------------------#
                           par(par.user)
                           leg.ncol    = min(3,pretty.box(n.now)$ncol)
                           leg.title   = expression(bold("Simulation"))
                           leg.legend  = if (any(now$parse)){
                                            parse(text=now$legend)
                                         }else{
                                            now$legend
                                         }#end if
                           xyz.plot( x              = x.list
                                   , y              = y.list
                                   , z              = z.list
                                   , fixed.xlim     = TRUE
                                   , fixed.ylim     = TRUE
                                   , xy.log         = xy.plog
                                   , pch            = pch.list
                                   , cex            = 1.2
                                   , nlevels        = n.colourbar
                                   , colour.palette = z.cscheme
                                   , xyz.title      = list( main      = letitre
                                                          , xlab      = lex
                                                          , ylab      = ley
                                                          , line.xlab = 4.1
                                                          , line.ylab = 2.6
                                                          , cex.main  = 1.1 * cex.ptsz
                                                          )#end list
                                   , xyz.sub        = pft.desc[pft.mp]
                                   , xyz.more       = list(grid=list(col=grid.colour
                                                                    ,lty="solid"))
                                   , key.log        = z.plog
                                   , key.title      = list(main=lacle,cex.main=.8*cex.main)
                                   , xyz.legend     = list( x       = "bottom"
                                                          , inset   = 0.0
                                                          , legend  = leg.legend
                                                          , pch     = now$pch
                                                          , col     = grey.fg
                                                          , border  = foreground
                                                          , bg      = background
                                                          , ncol    = leg.ncol
                                                          , title   = leg.title
                                                          , cex     = cex.ptsz
                                                          , xpd     = TRUE
                                                          )#end legend
                                   )#end xyz.plot
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           clean.tmp()
                           #---------------------------------------------------------------#
                        }#end for
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#











                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #     XYZ plots by DBH (annual means).                             #
                        #------------------------------------------------------------------#

                        #----- Split the arrays into lists. -------------------------------#
                        x.list    = split(x=x.tsdbh  ,f=ai.tsdbh[,4])
                        y.list    = split(x=y.tsdbh  ,f=ai.tsdbh[,4])
                        z.list    = split(x=z.tsdbh  ,f=ai.tsdbh[,4])
                        pch.list  = split(x=pch.tsdbh,f=ai.tsdbh[,4])
                        #------------------------------------------------------------------#



                        #------ Fix the titles. -------------------------------------------#
                        letitre   = paste("Annual means - ",z.desc,"\n",out.desc,sep="")
                        lex       = desc.unit(desc=x.desc,unit=x.unit)
                        ley       = desc.unit(desc=y.desc,unit=y.unit)
                        lacle     = desc.unit(desc=NULL  ,unit=z.unit)
                        #------------------------------------------------------------------#




                        #----- Plot the parameter space by DBH class. ---------------------#
                        for (o in 1:nout){
                           #----- Get the path. -------------------------------------------#
                           now.outpath  = outpath[[o]]$global[[g]]$scenario[[s]]$xyz.dbh[y]
                           #---------------------------------------------------------------#

                           #----- Open file or display. -----------------------------------#
                           fichier = paste(now.outpath,"/y-",y.vname,"_x-",x.vname
                                          ,"_z-",z.vname,"-",out.suffix,"-dbh."
                                          ,outform[o],sep="")
                           if (outform[o] == "x11"){
                              X11(width=size$width,height=size$height,pointsize=ptsz)
                           }else if(outform[o] == "png"){
                              png(filename=fichier,width=size$width*depth
                                 ,height=size$height*depth,pointsize=ptsz,res=depth)
                           }else if(outform[o] == "eps"){
                              postscript(file=fichier,width=size$width,height=size$height
                                        ,pointsize=ptsz,paper=size$paper)
                           }else if(outform[o] == "pdf"){
                              pdf(file=fichier,onefile=FALSE,width=size$width
                                 ,height=size$height,pointsize=ptsz,paper=size$paper)
                           }#end if
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #      Plot the boxes.                                          #
                           #---------------------------------------------------------------#
                           par(par.user)
                           leg.ncol    = min(3,pretty.box(n.now)$ncol)
                           leg.title   = expression(bold("Simulation"))
                           leg.legend  = if (any(now$parse)){
                                            parse(text=now$legend)
                                         }else{
                                            now$legend
                                         }#end if
                           xyz.plot( x              = x.list
                                   , y              = y.list
                                   , z              = z.list
                                   , fixed.xlim     = TRUE
                                   , fixed.ylim     = TRUE
                                   , xy.log         = xy.plog.dbh
                                   , pch            = pch.list
                                   , cex            = 1.2
                                   , nlevels        = n.colourbar
                                   , colour.palette = z.cscheme
                                   , xyz.title      = list( main      = letitre
                                                          , xlab      = lex
                                                          , ylab      = ley
                                                          , line.xlab = 4.1
                                                          , line.ylab = 2.6
                                                          , cex.main  = 1.1 * cex.ptsz
                                                          )#end list
                                   , xyz.sub        = dbh.desc[dbh.mp]
                                   , xyz.more       = list(grid=list(col=grid.colour
                                                                    ,lty="solid"))
                                   , key.log        = z.plog.dbh
                                   , key.title      = list(main=lacle,cex.main=.8*cex.main)
                                   , xyz.legend     = list( x       = "bottom"
                                                          , inset   = 0.0
                                                          , legend  = leg.legend
                                                          , pch     = now$pch
                                                          , col     = grey.fg
                                                          , border  = foreground
                                                          , bg      = background
                                                          , ncol    = leg.ncol
                                                          , title   = leg.title
                                                          , cex     = cex.ptsz
                                                          , xpd     = TRUE
                                                          )#end legend
                                   )#end xyz.plot
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           clean.tmp()
                           #---------------------------------------------------------------#
                        }#end for
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#














                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #     XYZ plots by PFT (DBH class, seasonal means).                #
                        #------------------------------------------------------------------#
                        for (e in 1:n.season){
                           for (d in 1:n.dbh){

                              #----- Keep only the elements that will go to this plot. ----#
                              x.tsnow = x.tsall[,,,e,d,]
                              y.tsnow = y.tsall[,,,e,d,]
                              z.tsnow = z.tsall[,,,e,d,]

                              #----- Create indices to split the arrays into lists. -------#
                              ai.tsnow = arrayInd(sequence(length(y.tsnow))
                                                 ,.dim=dim(y.tsnow))
                              #------------------------------------------------------------#


                              #----- Find the format of the points for the plots. ---------#
                              idx       = arrayInd(1:length(y.tsnow),.dim=dim(y.tsnow))[,1]
                              pch.tsnow = array(data=now$pch[idx],dim=dim(y.tsnow))
                              #------------------------------------------------------------#



                              #----- Split the arrays into lists. -------------------------#
                              x.list    = split(x=x.tsnow  ,f=ai.tsnow[,4])
                              y.list    = split(x=y.tsnow  ,f=ai.tsnow[,4])
                              z.list    = split(x=z.tsnow  ,f=ai.tsnow[,4])
                              pch.list  = split(x=pch.tsnow,f=ai.tsnow[,4])
                              #------------------------------------------------------------#



                              #------ Fix the titles. -------------------------------------#
                              letitre   = paste(z.desc," - ",season.mp.desc[e]," means"
                                               ," - DBH: ",dbh.desc[d]
                                               ,"\n",out.desc,sep="")
                              lex       = desc.unit(desc=x.desc,unit=x.unit)
                              ley       = desc.unit(desc=y.desc,unit=y.unit)
                              lacle     = desc.unit(desc=NULL  ,unit=z.unit)
                              #------------------------------------------------------------#




                              #----- Plot the parameter space by DBH class. ---------------#
                              for (o in 1:nout){
                                 #----- Get the path. -------------------------------------#
                                 scen.out     = outpath[[o]]$global[[g]]$scenario[[s]]
                                 now.outpath  = scen.out$xyz.all[d,e,y]
                                 #---------------------------------------------------------#

                                 #----- Open file or display. -----------------------------#
                                 fichier = paste(now.outpath,"/y-",y.vname,"_x-",x.vname
                                                ,"_z-",z.vname,"-",out.suffix,"-dbh_"
                                                ,dbh.key[d],"-season_",season.suffix[e],"."
                                                ,outform[o],sep="")
                                 if (outform[o] == "x11"){
                                    X11(width=size$width,height=size$height,pointsize=ptsz)
                                 }else if(outform[o] == "png"){
                                    png(filename=fichier,width=size$width*depth
                                       ,height=size$height*depth,pointsize=ptsz,res=depth)
                                 }else if(outform[o] == "eps"){
                                    postscript(file=fichier,width=size$width
                                              ,height=size$height,pointsize=ptsz
                                              ,paper=size$paper)
                                 }else if(outform[o] == "pdf"){
                                    pdf(file=fichier,onefile=FALSE,width=size$width
                                       ,height=size$height,pointsize=ptsz,paper=size$paper)
                                 }#end if
                                 #---------------------------------------------------------#



                                 #---------------------------------------------------------#
                                 #      Plot the boxes.                                    #
                                 #---------------------------------------------------------#
                                 par(par.user)
                                 leg.ncol    = min(3,pretty.box(n.now)$ncol)
                                 leg.title   = expression(bold("Simulation"))
                                 leg.legend  = if (any(now$parse)){
                                                  parse(text=now$legend)
                                               }else{
                                                  now$legend
                                               }#end if
                                 xyz.plot( x              = x.list
                                         , y              = y.list
                                         , z              = z.list
                                         , fixed.xlim     = TRUE
                                         , fixed.ylim     = TRUE
                                         , xy.log         = xy.plog
                                         , pch            = pch.list
                                         , cex            = 1.2
                                         , nlevels        = n.colourbar
                                         , colour.palette = z.cscheme
                                         , xyz.title      = list( main      = letitre
                                                                , xlab      = lex
                                                                , ylab      = ley
                                                                , line.xlab = 4.1
                                                                , line.ylab = 2.6
                                                                , cex.main  = 0.9*cex.ptsz
                                                                )#end list
                                         , xyz.sub        = pft.desc[pft.mp]
                                         , xyz.more       = list(grid=list(col=grid.colour
                                                                          ,lty="solid"))
                                         , key.log        = z.plog.dbh
                                         , key.title      = list(main=lacle
                                                                ,cex.main=0.8*cex.main)
                                         , xyz.legend     = list( x       = "bottom"
                                                                , inset   = 0.0
                                                                , legend  = leg.legend
                                                                , pch     = now$pch
                                                                , col     = grey.fg
                                                                , border  = foreground
                                                                , bg      = background
                                                                , ncol    = leg.ncol
                                                                , title   = leg.title
                                                                , cex     = cex.ptsz
                                                                , xpd     = TRUE
                                                                )#end legend
                                         )#end xyz.plot
                                 #---------------------------------------------------------#



                                 #----- Close the device. ---------------------------------#
                                 if (outform[o] == "x11"){
                                    locator(n=1)
                                    dev.off()
                                 }else{
                                    dev.off()
                                 }#end if
                                 clean.tmp()
                                 #---------------------------------------------------------#
                              }#end for (o in 1:nout)
                              #------------------------------------------------------------#
                           }#end for (d in 1:n.dbh)
                           #---------------------------------------------------------------#
                        }#end for (e in 1:n.season.mp)
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                     }#end if (! z.vname in c(x.vname,y.vname))
                     #---------------------------------------------------------------------#
                  }#end for (z in loop.z)
                  #------------------------------------------------------------------------#
               }#end if (x.vname != y.vname)
               #---------------------------------------------------------------------------#
            }#end for (x in loop.x)
            #------------------------------------------------------------------------------#
         }#end for (y in loop.y)
         #.................................................................................#
         #.................................................................................#
      }#end for (s in loop.scenario)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   }#end for (p in loop.panel)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#




   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #        Now we plot panels comparing the different scenarios.  We use the same list    #
   # of variables as the time series.                                                      #
   #---------------------------------------------------------------------------------------#
   if (n.scenario == 2){
      if (is.na(plot.scencomp)){
         loop.scencomp = sequence(nscen.comp)[1]
      }else if (plot.scencomp){
         loop.scencomp = sequence(nscen.comp)
      }else{
         loop.scencomp = integer(0)
      }#end if
   }else{
      loop.scencomp = integer(0)
   }#end if
   cat("     * Plot the scenario comparison...","\n")
   for (v in loop.scencomp){
      #----- Copy variable info. ----------------------------------------------------------#
      var.vname   = scen.comp$vname  [v]
      var.desc    = scen.comp$desc   [v]
      var.quant   = scen.comp$quant  [v]
      var.low     = scen.comp$low    [v]
      var.high    = scen.comp$high   [v]
      var.unit    = scen.comp$unit   [v]
      is.pft      = scen.comp$pft    [v]
      is.dbh      = scen.comp$dbh    [v]
      var.plt     = scen.comp$plt    [v]
      var.pft     = paste(var.vname,"pft"   ,sep="")
      var.pftdbh  = paste(var.vname,"pftdbh",sep="")
      #------------------------------------------------------------------------------------#



      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
      #     Decide whether to plot or not.                                                 #
      #------------------------------------------------------------------------------------#
      if (var.plt){
         cat  ("       . ",var.desc,"...","\n")
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the variable difference by panel.                                     #
         #---------------------------------------------------------------------------------#



         #----- Retrieve all simulations from this panel. ---------------------------------#
         this    = eft[[var.vname]]$ts[,,y.sel,,drop=FALSE]
         this    = apply(X = this, FUN = mean, MARGIN = c(1,4),na.rm=TRUE)
         o.dim   = dim(this)[-1]      ; names(o.dim  ) = "season"
         o.dname = dimnames(this)[-1] ; names(o.dname) = "season"
         o.seq   = mapply( FUN = sequence, nvec = o.dim, SIMPLIFY = FALSE)
         #---------------------------------------------------------------------------------#



         #----- Find the panel information. -----------------------------------------------#
         p.type         = simul$dim.type == "panel"
         p.name         = "panel"
         p.ind          = simul$panel$index  [g.sel]
         p.key          = simul$key[g.sel,p.type]
         p.level        = list(simul$panel$level)
         p.dim          = simul$panel$n.level
         p.seq          = sequence(p.dim)
         names(p.level) = "panel"
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Retrieve the scenarios and the default for each scenario.                  #
         #---------------------------------------------------------------------------------#
         s.type         = simul$dim.type == "scenario"
         s.name         = names(simul$dim[s.type])
         s.dim          = simul$dim[s.type]
         s.dnames       = simul$panel$level
         s.value        = simul$value[g.sel,s.type]
         s.key          = simul$key[g.sel,s.type]
         s.level        = split(x=simul$scenario$level,f=simul$scenario$idxcol)
         names(s.level) = names(simul$dim)[s.type]
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the array names, and determine the index that has the "zero" scenario #
         # so the figures are relative to the reference scenario.                          #
         #---------------------------------------------------------------------------------#
         if (n.scenario == 1){
            s.ind           = simul$scenario$index[g.sel]
            zero            = s.key == simul$default[s.type]
            a.dim           = sapply(X=c(s.dim,o.dim),FUN=c)
            a.dnames        = c(list(s.level),o.dname)
            names(a.dnames) = c(names(simul$dim)[s.type],names(o.dim))
         }else{
            s.default = split(simul$default[s.type],index(simul$default[s.type]))
            s.ind     = data.frame( apply(X=simul$scenario$index,MARGIN=2,FUN="[",g.sel)
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
            zero      = apply( X      = mapply(FUN="==", s.key,s.default)
                             , MARGIN = 1
                             , FUN    = all
                             )#end apply
            a.dim     = sapply(X=c(p.dim,s.dim  ,o.dim  ),FUN=c)
            a.dnames  = c(p.level,s.level,o.dname)
            x.zero    = scenario[[1]]$value[match(s.default[[1]],scenario[[1]]$key)]
            y.zero    = scenario[[2]]$value[match(s.default[[2]],scenario[[2]]$key)]
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Build the lists to split the array by simulation and panel.                #
         #---------------------------------------------------------------------------------#
         if (! is.list(p.ind)) p.ind = list(p.ind)
         if (! is.list(s.ind)) s.ind = list(s.ind)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the zero for all simulations.                                         #
         #---------------------------------------------------------------------------------#
         zero = as.matrix(expand.grid(c(list(rep(which(zero), times=prod(s.dim))),o.seq)))
         this = this - this[zero]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Split the stuff into panels, and create the scenario maps.                 #
         #---------------------------------------------------------------------------------#
         this       = array(data=this,dim=a.dim,dimnames=a.dnames)
         n.a.dim    = length(a.dim)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Plot the data by season.                                                   #
         #---------------------------------------------------------------------------------#
         for (e in 1:n.season){

            #------------------------------------------------------------------------------#
            #     Keep the current season only.                                            #
            #------------------------------------------------------------------------------#
            bye      = n.a.dim
            inds     = as.matrix(expand.grid(mapply(FUN=sequence,a.dim[-bye])))
            inds     = cbind(inds,e)
            this.var = array(this[inds],dim=a.dim[-bye],dimnames=a.dnames[-bye])
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Set the annotation and some auxiliary variables.  These depend on the   #
            # number of scenario dimensions.                                               #
            #------------------------------------------------------------------------------#
            letitre = paste("Scenario comparision - ",var.desc,"\n",season.desc[e]
                           ," - ",global.desc,sep="")
            if (n.scenario == 1){
               now    = scenario[[1]]
               n.now  = length(now$key)


               #---------------------------------------------------------------------------#
               #      Find the axis limits and labels.                                     #
               #---------------------------------------------------------------------------#
               xlimit  = pretty.xylim(u=now$value,fracexp=0.,is.log=FALSE)
               ylimit  = pretty.xylim(u=this.var ,fracexp=0.,is.log=FALSE)
               x.value = now$value
               x.label = now$label
               lex     = now$alabel
               ley     = desc.unit(desc=var.desc,unit=var.unit)
               #---------------------------------------------------------------------------#

            }else{
               #----- Axis labels. --------------------------------------------------------#
               lex     = scenario[[1]]$alabel
               x.value = scenario[[1]]$value
               x.label = scenario[[1]]$label
               ley     = scenario[[2]]$alabel
               y.value = scenario[[2]]$value
               y.label = scenario[[2]]$label
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    Break the values into bins.                                            #
               #---------------------------------------------------------------------------#
               var.scheme = two.palettes( x     = this.var
                                        , n     = n.colourbar
                                        , white = n.whitebar
                                        , low   = var.low
                                        , high  = var.high
                                        )#end two.palettes
               var.brks    = var.scheme$breaks
               var.colours = var.scheme$colours
               n.brks      = var.scheme$n.breaks
               #---------------------------------------------------------------------------#


               #----- Make the edges. -----------------------------------------------------#
               xleft       = var.brks[-n.brks]
               xright      = var.brks[     -1]
               ybottom     = rep(0,times=n.brks)
               ytop        = rep(1,times=n.brks)
               legat       = pretty(var.brks,n=10)
               legat       = legat[legat >= min(var.brks) & legat <= max(var.brks)]
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Loop over all output formats.                                            #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Get the path and create file name. ----------------------------------#
               if ( e == n.season){
                  now.outpath  = outpath[[o]]$global[[g]]$scenpanel.year
                  fichier      = paste(now.outpath,"/scencomp-year-"
                                      ,var.vname,"-",global.suffix,".",outform[o],sep="")
               }else{
                  now.outpath  = outpath[[o]]$global[[g]]$scenpanel.season[e]
                  fichier      = paste(now.outpath,"/scencomp-",season.suffix[e],"-"
                                      ,var.vname,"-",global.suffix,".",outform[o],sep="")
               }#end if
               #---------------------------------------------------------------------------#


               #----- Open the file. ------------------------------------------------------#
               if (outform[o] == "x11"){
                  X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=wsize$width*depth
                     ,height=wsize$height*depth,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=wsize$width,height=wsize$height
                            ,pointsize=ptsz,paper=wsize$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=wsize$width
                     ,height=wsize$height,pointsize=ptsz,paper=wsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Save the margins to avoid losing the data. --------------------------#
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Split the plotting window.                                           #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma = c(0.2,5.1,4.5,0))
               layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
                     , heights = c(rep(9/lo.panel$nrow,lo.panel$nrow),2)
                     )#end layout
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Plot the legend.  What goes into the legend depends on the number of  #
               # scenarios.                                                                #
               #---------------------------------------------------------------------------#
               if (n.scenario == 1){
                  par(mar=c(0.2,4.4,0.1,2.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1))
                  legend( x      = "bottom"
                        , inset  = inset
                        , legend = if (any(now$parse)){
                                      parse(text=now$legend)
                                   }else{
                                      now$legend
                                   }#end if
                        , col    = now$colour
                        , pch    = 16
                        , lwd    = 2.5
                        , title  = expression(bold(now$alabel))
                        , cex    = 1.0
                        , xpd    = TRUE
                        )#end legend
               }else{
                  par(mar=c(3,3,2,3)+0.1)
                  plot.new()
                  plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop)
                             ,xaxs="i",yaxs="i")
                  rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=var.colours)
                  box()
                  axis.rt(side=1,at=legat,las=5,off=0.5)
                  title(main=desc.unit(desc=paste(var.desc,"change"),unit=var.unit)
                       ,xlab="",ylab="",cex.main=cex.ptsz)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all panels, and plot the data for this season.             #
               #---------------------------------------------------------------------------#
               loop.panel = sequence(max(1,simul$panel$n.level   ))
               for (p in loop.panel){

                  #----- Find out where is this box going, and set up axes and margins. ---#
                  left    = (p %% lo.panel$ncol) == 1 || lo.panel$ncol == 1
                  right   = (p %% lo.panel$ncol) == 0
                  top     = p <= lo.panel$ncol
                  bottom  = p > (lo.panel$nrow - 1) * lo.panel$ncol
                  if (n.scenario == 1){
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                  }else{
                     mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Make the sub-title. ----------------------------------------------#
                  lesub=simul$panel$title[p]
                  #------------------------------------------------------------------------#



                  #----- Set the window. --------------------------------------------------#
                  par(mar = mar.now)
                  if (n.scenario == 1){
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit)
                     if (bottom) axis(side=1,at=x.value,labels=x.label)
                     if (left)   axis(side=2)
                     box()
                     grid(col=grid.colour,lty="solid")
                     title(main=lesub,cex.main=cex.ptsz,line=0.5)
                     abline(h=0,v=0,col=foreground,lty="dotdash",lwd=2.0)
                     for (n in 1:n.now){
                        points(x = now$value,y=this.var[p,],type="o",col=now$colour[n]
                              ,pch=16,lwd=2.5)
                     }#end for
                  }else{
                     image(x=x.value,y=y.value,z=this.var[p,,],col=var.colours
                          ,breaks=var.brks
                          ,xaxt="n",yaxt="n",main=lesub,cex.main=cex.ptsz,xlab="",ylab="")
                     text (x=x.zero,y=y.zero,labels="0",col=foreground,font=2,cex=2)
                     if (bottom) axis(side=1,at=x.value,labels=x.label,cex.axis=1.1,las=1)
                     if (left)   axis(side=2,at=y.value,labels=y.label,cex.axis=1.1,las=1)
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               gtitle( main      = letitre
                     , xlab      = lex
                     , ylab      = ley
                     , off.xlab  = 2/11
                     , line.xlab = 4.1
                     , line.ylab = 2.6
                     , cex.main  = 1.1*cex.ptsz
                     )#end gtitle
               # par(las=0)
               # if (n.scenario == 1){
               #    mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
               #    mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
               # }else{
               #    mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff.im)
               #    mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff.im)
               # }#end if
               # mtext(side=3,text=letitre,outer=TRUE,cex=1.1*cex.ptsz,font=2)
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               clean.tmp()
               #---------------------------------------------------------------------------#
            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end for (e in 1:n.season)
         #---------------------------------------------------------------------------------#
      }#end if
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#






      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
      #     Check whether to plot the comparison by PFT.                                   #
      #------------------------------------------------------------------------------------#
      if (var.plt && is.pft){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the variable difference by panel.                                     #
         #---------------------------------------------------------------------------------#



         #----- Retrieve all simulations from this panel. ---------------------------------#
         this    = eft[[var.vname]]$tspft[,,y.sel,,]
         this    = apply(X = this, FUN = mean, MARGIN = c(1,4,5),na.rm=TRUE)
         o.dim   = dim(this)[-1]      ; names(o.dim  ) = c("season","pft")
         o.dname = dimnames(this)[-1] ; names(o.dname) = c("season","pft")
         o.seq   = mapply( FUN = sequence, nvec = o.dim, SIMPLIFY = FALSE)
         #---------------------------------------------------------------------------------#



         #----- Find the panel information. -----------------------------------------------#
         p.type         = simul$dim.type == "panel"
         p.name         = "panel"
         p.ind          = simul$panel$index  [g.sel]
         p.key          = simul$key[g.sel,p.type]
         p.level        = list(simul$panel$level)
         p.dim          = simul$panel$n.level
         p.seq          = sequence(p.dim)
         names(p.level) = "panel"
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Retrieve the scenarios and the default for each scenario.                  #
         #---------------------------------------------------------------------------------#
         s.type         = simul$dim.type == "scenario"
         s.name         = names(simul$dim[s.type])
         s.dim          = simul$dim[s.type]
         s.dnames       = simul$panel$level
         s.value        = simul$value[g.sel,s.type]
         s.key          = simul$key[g.sel,s.type]
         s.level        = split(x=simul$scenario$level,f=simul$scenario$idxcol)
         names(s.level) = names(simul$dim)[s.type]
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the array names, and determine the index that has the "zero" scenario #
         # so the figures are relative to the reference scenario.                          #
         #---------------------------------------------------------------------------------#
         if (n.scenario == 1){
            s.ind           = simul$scenario$index[g.sel]
            zero            = s.key == simul$default[s.type]
            a.dim           = sapply(X=c(s.dim,o.dim),FUN=c)
            a.dnames        = c(list(s.level),o.dname)
            names(a.dnames) = c(names(simul$dim)[s.type],names(o.dim))
         }else{
            s.default = split(simul$default[s.type],index(simul$default[s.type]))
            s.ind     = data.frame( apply(X=simul$scenario$index,MARGIN=2,FUN="[",g.sel)
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
            zero      = apply( X      = mapply(FUN="==", s.key,s.default)
                             , MARGIN = 1
                             , FUN    = all
                             )#end apply
            a.dim     = sapply(X=c(p.dim,s.dim  ,o.dim  ),FUN=c)
            a.dnames  = c(p.level,s.level,o.dname)
            x.zero    = scenario[[1]]$value[match(s.default[[1]],scenario[[1]]$key)]
            y.zero    = scenario[[2]]$value[match(s.default[[2]],scenario[[2]]$key)]
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Build the lists to split the array by simulation and panel.                #
         #---------------------------------------------------------------------------------#
         if (! is.list(p.ind)) p.ind = list(p.ind)
         if (! is.list(s.ind)) s.ind = list(s.ind)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the zero for all simulations.                                         #
         #---------------------------------------------------------------------------------#
         zero = as.matrix(expand.grid(c(list(rep(which(zero), times=prod(s.dim))),o.seq)))
         this = this - this[zero]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Split the stuff into panels, and create the scenario maps.                 #
         #---------------------------------------------------------------------------------#
         this       = array(data=this,dim=a.dim,dimnames=a.dnames)
         n.a.dim    = length(a.dim)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Plot the data by PFT and season.                                           #
         #---------------------------------------------------------------------------------#
         for (f in 1:(n.pft-1)){
            for (e in 1:n.season){

               #------ Keep only the current season. --------------------------------------#
               bye      = n.a.dim + c(-1,0)
               inds     = as.matrix(expand.grid(mapply(FUN=sequence,a.dim[-bye])))
               inds     = cbind(inds,e,f)
               this.var = array(this[inds],dim=a.dim[-bye],dimnames=a.dnames[-bye])
               #---------------------------------------------------------------------------#


               
               #---------------------------------------------------------------------------#
               #      Set the annotation and some auxiliary variables.  These depend on    #
               # the number of scenario dimensions.                                        #
               #---------------------------------------------------------------------------#
               letitre = paste("Scenario comparision - ",var.desc
                              ,"\n",pft.desc[f]," - ",season.desc[e]
                              ,"\n",global.desc,sep="")
               if (n.scenario == 1){
                  now    = scenario[[1]]
                  n.now  = length(now$key)


                  #------------------------------------------------------------------------#
                  #      Find the axis limits and labels.                                  #
                  #------------------------------------------------------------------------#
                  xlimit  = pretty.xylim(u=now$value,fracexp=0.,is.log=FALSE)
                  ylimit  = pretty.xylim(u=this.var,fracexp=0.,is.log=FALSE)
                  lex     = now$alabel
                  ley     = desc.unit(desc=var.desc,unit=var.unit)
                  x.value = now$value
                  x.label = now$label
                  #------------------------------------------------------------------------#

               }else{
                  #----- Axis labels. -----------------------------------------------------#
                  lex     = scenario[[1]]$alabel
                  x.value = scenario[[1]]$value
                  x.label = scenario[[1]]$label
                  ley     = scenario[[2]]$alabel
                  y.value = scenario[[2]]$value
                  y.label = scenario[[2]]$label
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #    Break the values into bins.                                         #
                  #------------------------------------------------------------------------#
                  var.scheme = two.palettes( x     = this.var
                                           , n     = n.colourbar
                                           , white = n.whitebar
                                           , low   = var.low
                                           , high  = var.high
                                           )#end two.palettes
                  var.brks    = var.scheme$breaks
                  var.colours = var.scheme$colours
                  n.brks      = var.scheme$n.breaks
                  #------------------------------------------------------------------------#


                  #----- Make the edges. --------------------------------------------------#
                  xleft       = var.brks[-n.brks]
                  xright      = var.brks[     -1]
                  ybottom     = rep(0,times=n.brks)
                  ytop        = rep(1,times=n.brks)
                  legat       = var.brks
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Loop over all output formats.                                         #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path and create file name. -------------------------------#
                  if ( e == n.season){
                     now.outpath  = outpath[[o]]$global[[g]]$scenpanelpft.year[f]
                     fichier      = paste(now.outpath,"/scencomp-year-",pft.suffix[f],"-"
                                         ,var.vname,"-",global.suffix,".",outform[o]
                                         ,sep="")
                  }else{
                     now.outpath  = outpath[[o]]$global[[g]]$scenpanelpft.season[f,e]
                     fichier      = paste(now.outpath,"/scencomp-",pft.suffix[f],"-"
                                         ,season.suffix[e],"-",var.vname,"-",global.suffix
                                         ,".",outform[o],sep="")
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Open the file. ---------------------------------------------------#
                  if (outform[o] == "x11"){
                     X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wsize$width*depth
                        ,height=wsize$height*depth,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wsize$width,height=wsize$height
                               ,pointsize=ptsz,paper=wsize$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=wsize$width
                        ,height=wsize$height,pointsize=ptsz,paper=wsize$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Save the margins to avoid losing the data. -----------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Split the plotting window.                                        #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par(oma = c(0.2,3,4.5,0))
                  layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
                        , heights = c(rep(9/lo.panel$nrow,lo.panel$nrow),2)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the legend.  What goes into the legend depends on the number  #
                  # of scenarios.                                                          #
                  #------------------------------------------------------------------------#
                  if (n.scenario == 1){
                     par(mar=c(0.2,4.4,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = inset
                           , legend = if (any(now$parse)){
                                         parse(text=now$legend)
                                      }else{
                                         now$legend
                                      }#end if
                           , col    = now$colour
                           , pch    = 16
                           , lwd    = 2.5
                           , title  = expression(bold(now$alabel))
                           , cex    = 1.0
                           , xpd    = TRUE
                           )#end legend
                  }else{
                     par(mar=c(3,3,2,3)+0.1)
                     plot.new()
                     plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop)
                                ,xaxs="i",yaxs="i")
                     rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop
                         ,col=var.colours)
                     box()
                     axis.rt(side=1,at=legat,las=5,off=0.5)
                     title(main=desc.unit(desc=paste(var.desc,"change"),unit=var.unit)
                          ,xlab="",ylab="",cex.main=cex.ptsz)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all panels, and plot the data for this season.          #
                  #------------------------------------------------------------------------#
                  loop.panel = sequence(max(1,simul$panel$n.level   ))
                  for (p in loop.panel){

                     #----- Find out where this box goes, and set up axes and margins. ----#
                     left    = (p %% lo.panel$ncol) == 1 || lo.panel$ncol == 1
                     right   = (p %% lo.panel$ncol) == 0
                     top     = p <= lo.panel$ncol
                     bottom  = p > (lo.panel$nrow - 1) * lo.panel$ncol
                     if (n.scenario == 1){
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     }else{
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     }#end if
                     #---------------------------------------------------------------------#


                     #----- Make the sub-title. -------------------------------------------#
                     lesub=simul$panel$title[p]
                     #---------------------------------------------------------------------#



                     #----- Set the window. -----------------------------------------------#
                     par(mar = mar.now)
                     if (n.scenario == 1){
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit)
                        if (bottom) axis(side=1,at=x.value,labels=x.label)
                        if (left)   axis(side=2)
                        box()
                        grid(col=grid.colour,lty="solid")
                        title(main=lesub,cex.main=cex.ptsz,line=0.5)
                        abline(h=0,v=0,col=foreground,lty="dotdash",lwd=2.0)
                        for (n in 1:n.now){
                           points(x = now$value,y=this.var[p,],type="o",col=now$colour[n]
                                 ,pch=16,lwd=2.5)
                        }#end for
                     }else{
                        image(x=scenario[[1]]$value,y=scenario[[2]]$value
                             ,z=this.var[p,,],col=var.colours,breaks=var.brks
                             ,xaxt="n",yaxt="n",main=lesub,xlab="",ylab=""
                             ,cex.main=cex.ptsz)
                        text (x=x.zero,y=y.zero,labels="0",col=foreground,font=2,cex=2)
                        if (bottom) axis(side=1,at=x.value,labels=x.label,las=1)
                        if (left)   axis(side=2,at=y.value,labels=y.label,las=1)
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 2/11
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # if (n.scenario == 1){
                  #    mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  #    mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # }else{
                  #    mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff.im)
                  #    mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff.im)
                  # }#end if
                  # mtext(side=3,text=letitre,outer=TRUE,cex=1.0,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #---------------------------------------------------------------------------#
            }#end for (e in 1:n.season)
            #------------------------------------------------------------------------------#
         }#end for (f in 1:(n.pft-1))
         #---------------------------------------------------------------------------------#
      }#end if (var.plt && is.pft)
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#






      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
      #     Check whether to plot the comparison by DBH.                                   #
      #------------------------------------------------------------------------------------#
      if (var.plt && is.pft && is.dbh){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the variable difference by panel.                                     #
         #---------------------------------------------------------------------------------#



         #----- Retrieve all simulations from this panel. ---------------------------------#
         this    = eft[[var.vname]]$tspftdbh[,,y.sel,,,n.pft]
         this    = apply(X = this, FUN = mean, MARGIN = c(1,4,5),na.rm=TRUE)
         o.dim   = dim(this)[-1]      ; names(o.dim  ) = c("season","dbh")
         o.dname = dimnames(this)[-1] ; names(o.dname) = c("season","dbh")
         o.seq   = mapply( FUN = sequence, nvec = o.dim, SIMPLIFY = FALSE)
         #---------------------------------------------------------------------------------#



         #----- Find the panel information. -----------------------------------------------#
         p.type         = simul$dim.type == "panel"
         p.name         = "panel"
         p.ind          = simul$panel$index  [g.sel]
         p.key          = simul$key[g.sel,p.type]
         p.level        = list(simul$panel$level)
         p.dim          = simul$panel$n.level
         p.seq          = sequence(p.dim)
         names(p.level) = "panel"
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Retrieve the scenarios and the default for each scenario.                  #
         #---------------------------------------------------------------------------------#
         s.type         = simul$dim.type == "scenario"
         s.name         = names(simul$dim[s.type])
         s.dim          = simul$dim[s.type]
         s.dnames       = simul$panel$level
         s.value        = simul$value[g.sel,s.type]
         s.key          = simul$key[g.sel,s.type]
         s.level        = split(x=simul$scenario$level,f=simul$scenario$idxcol)
         names(s.level) = names(simul$dim)[s.type]
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the array names, and determine the index that has the "zero" scenario #
         # so the figures are relative to the reference scenario.                          #
         #---------------------------------------------------------------------------------#
         if (n.scenario == 1){
            s.ind           = simul$scenario$index[g.sel]
            zero            = s.key == simul$default[s.type]
            a.dim           = sapply(X=c(s.dim,o.dim),FUN=c)
            a.dnames        = c(list(s.level),o.dname)
            names(a.dnames) = c(names(simul$dim)[s.type],names(o.dim))
         }else{
            s.default = split(simul$default[s.type],index(simul$default[s.type]))
            s.ind     = data.frame( apply(X=simul$scenario$index,MARGIN=2,FUN="[",g.sel)
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
            zero      = apply( X      = mapply(FUN="==", s.key,s.default)
                             , MARGIN = 1
                             , FUN    = all
                             )#end apply
            a.dim     = sapply(X=c(p.dim,s.dim  ,o.dim  ),FUN=c)
            a.dnames  = c(p.level,s.level,o.dname)
            x.zero    = scenario[[1]]$value[match(s.default[[1]],scenario[[1]]$key)]
            y.zero    = scenario[[2]]$value[match(s.default[[2]],scenario[[2]]$key)]
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Build the lists to split the array by simulation and panel.                #
         #---------------------------------------------------------------------------------#
         if (! is.list(p.ind)) p.ind = list(p.ind)
         if (! is.list(s.ind)) s.ind = list(s.ind)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the zero for all simulations.                                         #
         #---------------------------------------------------------------------------------#
         zero = as.matrix(expand.grid(c(list(rep(which(zero), times=prod(s.dim))),o.seq)))
         this = this - this[zero]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Split the stuff into panels, and create the scenario maps.                 #
         #---------------------------------------------------------------------------------#
         this       = array(data=this,dim=a.dim,dimnames=a.dnames)
         n.a.dim    = length(a.dim)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Plot the data by PFT and season.                                           #
         #---------------------------------------------------------------------------------#
         for (d in 1:n.dbh){
            for (e in 1:n.season){

               #------ Keep only the current season. --------------------------------------#
               bye      = n.a.dim + c(-1,0)
               inds     = as.matrix(expand.grid(mapply(FUN=sequence,a.dim[-bye])))
               inds     = cbind(inds,e,f)
               this.var = array(this[inds],dim=a.dim[-bye],dimnames=a.dnames[-bye])
               #---------------------------------------------------------------------------#


               
               #---------------------------------------------------------------------------#
               #      Set the annotation and some auxiliary variables.  These depend on    #
               # the number of scenario dimensions.                                        #
               #---------------------------------------------------------------------------#
               letitre = paste("Change in mean ",var.desc
                              ,"\n","DBH class: ",dbh.desc[d]," - ",season.desc[e]
                              ,"\n",global.desc,sep="")
               if (n.scenario == 1){
                  now    = scenario[[1]]
                  n.now  = length(now$key)


                  #------------------------------------------------------------------------#
                  #      Find the axis limits and labels.                                  #
                  #------------------------------------------------------------------------#
                  xlimit  = pretty.xylim(u=now$value,fracexp=0.,is.log=FALSE)
                  ylimit  = pretty.xylim(u=this.var,fracexp=0.,is.log=FALSE)
                  x.value = now$value
                  x.label = now$label
                  lex     = now$alabel
                  ley     = desc.unit(desc=var.desc,unit=var.unit)
                  #------------------------------------------------------------------------#

               }else{
                  #----- Axis labels. -----------------------------------------------------#
                  lex     = scenario[[1]]$alabel
                  x.value = scenario[[1]]$value
                  x.label = scenario[[1]]$label
                  ley     = scenario[[2]]$alabel
                  y.value = scenario[[2]]$value
                  y.label = scenario[[2]]$label
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #    Break the values into bins.                                         #
                  #------------------------------------------------------------------------#
                  var.scheme = two.palettes( x     = this.var
                                           , n     = n.colourbar
                                           , white = n.whitebar
                                           , low   = var.low
                                           , high  = var.high
                                           )#end two.palettes
                  var.brks    = var.scheme$breaks
                  var.colours = var.scheme$colours
                  n.brks      = var.scheme$n.breaks
                  #------------------------------------------------------------------------#


                  #----- Make the edges. --------------------------------------------------#
                  xleft       = var.brks[-n.brks]
                  xright      = var.brks[     -1]
                  ybottom     = rep(0,times=n.brks)
                  ytop        = rep(1,times=n.brks)
                  legat       = var.brks
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Loop over all output formats.                                         #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Get the path and create file name. -------------------------------#
                  if ( e == n.season){
                     now.outpath  = outpath[[o]]$global[[g]]$scenpaneldbh.year[d]
                     fichier      = paste(now.outpath,"/scencomp-year-",dbh.suffix[d],"-"
                                         ,var.vname,"-",global.suffix,".",outform[o]
                                         ,sep="")
                  }else{
                     now.outpath  = outpath[[o]]$global[[g]]$scenpaneldbh.season[d,e]
                     fichier      = paste(now.outpath,"/scencomp-",dbh.suffix[d],"-"
                                         ,season.suffix[e],"-",var.vname,"-",global.suffix
                                         ,".",outform[o],sep="")
                  }#end if
                  #------------------------------------------------------------------------#


                  #----- Open the file. ---------------------------------------------------#
                  if (outform[o] == "x11"){
                     X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wsize$width*depth
                        ,height=wsize$height*depth,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wsize$width,height=wsize$height
                               ,pointsize=ptsz,paper=wsize$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE,width=wsize$width
                        ,height=wsize$height,pointsize=ptsz,paper=wsize$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Save the margins to avoid losing the data. -----------------------#
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Split the plotting window.                                        #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par(oma = c(0.2,3,4.5,0))
                  layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
                        , heights = c(rep(9/lo.panel$nrow,lo.panel$nrow),2)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the legend.  What goes into the legend depends on the number  #
                  # of scenarios.                                                          #
                  #------------------------------------------------------------------------#
                  if (n.scenario == 1){
                     par(mar=c(0.2,4.4,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = inset
                           , legend = if (any(now$parse)){
                                         parse(text=now$legend)
                                      }else{
                                         now$legend
                                      }#end if
                           , col    = now$colour
                           , pch    = 16
                           , lwd    = 2.5
                           , title  = expression(bold(now$alabel))
                           , cex    = 1.0
                           , xpd    = TRUE
                           )#end legend
                  }else{
                     par(mar=c(3,3,2,3)+0.1)
                     plot.new()
                     plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop)
                                ,xaxs="i",yaxs="i")
                     rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop
                         ,col=var.colours)
                     box()
                     axis.rt(side=1,at=legat,las=5,off=0.5)
                     title(main=desc.unit(desc=paste(var.desc,"change"),unit=var.unit)
                          ,xlab="",ylab="",cex.main=cex.ptsz)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all panels, and plot the data for this season.          #
                  #------------------------------------------------------------------------#
                  loop.panel = sequence(max(1,simul$panel$n.level   ))
                  for (p in loop.panel){

                     #----- Find out where this box goes, and set up axes and margins. ----#
                     left    = (p %% lo.panel$ncol) == 1 || lo.panel$ncol == 1
                     right   = (p %% lo.panel$ncol) == 0
                     top     = p <= lo.panel$ncol
                     bottom  = p > (lo.panel$nrow - 1) * lo.panel$ncol
                     if (n.scenario == 1){
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     }else{
                        mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                     }#end if
                     #---------------------------------------------------------------------#


                     #----- Make the sub-title. -------------------------------------------#
                     lesub=simul$panel$title[p]
                     #---------------------------------------------------------------------#



                     #----- Set the window. -----------------------------------------------#
                     par(mar = mar.now)
                     if (n.scenario == 1){
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit)
                        if (bottom) axis(side=1,at=x.value,labels=x.label)
                        if (left)   axis(side=2)
                        box()
                        grid(col=grid.colour,lty="solid")
                        title(main=lesub,cex.main=cex.ptsz,line=0.5)
                        abline(h=0,v=0,col=foreground,lty="dotdash",lwd=2.0)
                        for (n in 1:n.now){
                           points(x = now$value,y=this.var[p,],type="o",col=now$colour[n]
                                 ,pch=16,lwd=2.5)
                        }#end forscen
                     }else{
                        image(x=scenario[[1]]$value,y=scenario[[2]]$value
                             ,z=this.var[p,,],col=var.colours,breaks=var.brks
                             ,xaxt="n",yaxt="n",main=lesub,xlab="",ylab=""
                             ,cex.main=cex.ptsz)
                        text (x=x.zero,y=y.zero,labels="0",col=foreground,font=2,cex=2)
                        if (bottom) axis(side=1,at=x.value,labels=x.label,las=1)
                        if (left)   axis(side=2,at=y.value,labels=y.label,las=1)
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  gtitle( main      = letitre
                        , xlab      = lex
                        , ylab      = ley
                        , off.xlab  = 2/11
                        , line.xlab = 4.1
                        , line.ylab = 2.6
                        , cex.main  = 1.1*cex.ptsz
                        )#end gtitle
                  # par(las=0)
                  # if (n.scenario == 1){
                  #    mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  #    mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  # }else{
                  #    mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff.im)
                  #    mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff.im)
                  # }#end if
                  # mtext(side=3,text=letitre,outer=TRUE,cex=1.0,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #---------------------------------------------------------------------------#
            }#end for (e in 1:n.season)
            #------------------------------------------------------------------------------#
         }#end for (d in 1:n.dbh)
         #---------------------------------------------------------------------------------#
      }#end if (var.plt && is.pft && is.dbh)
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
   }#end for (v in 1:nscen.comp)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#




   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #        Now we plot bo plots comparing the different scenarios.  We use the same list  #
   # of variables as the time series.                                                      #
   #---------------------------------------------------------------------------------------#
   if (n.scenario == 2){
      if (is.na(plot.panelbox)){
         loop.panelbox = sequence(npanel.box)[1]
      }else if (plot.panelbox){
         loop.panelbox = sequence(npanel.box)
      }else{
         loop.panelbox = integer(0)
      }#end if
   }else{
      loop.panelbox = integer(0)
   }#end if
   cat("     * Plot the multiple box plot...","\n")
   for (v in loop.panelbox){
      #----- Copy variable info. ----------------------------------------------------------#
      var.vname    = panel.box$vname  [v]
      var.desc     = panel.box$desc   [v]
      var.quant    = panel.box$quant  [v]
      var.low      = panel.box$low    [v]
      var.high     = panel.box$high   [v]
      var.unit     = panel.box$unit   [v]
      var.plog     = panel.box$plog   [v]
      is.pft       = panel.box$pft    [v]
      is.dbh       = panel.box$dbh    [v]
      var.plt      = panel.box$plt    [v]
      var.pft      = paste(var.vname,"pft"   ,sep="")
      var.pftdbh   = paste(var.vname,"pftdbh",sep="")
      if (var.plog){
         var.xylog = "y"
      }else{
         var.xylog = ""
      }#end if
      #------------------------------------------------------------------------------------#



      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
      #     Decide whether to plot or not.                                                 #
      #------------------------------------------------------------------------------------#
      if (var.plt){
         cat  ("       . ",var.desc,"...","\n")
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the box plot by panel.                                                #
         #---------------------------------------------------------------------------------#



         #----- Retrieve all simulations from this panel. ---------------------------------#
         if (boxplot.aggr %in% c("unlist","c")){
            original = eft[[var.vname]]$ts[,,y.sel,,drop=FALSE]
            concat   = apply( X        = expand.grid( dimnames(original)[[2]]
                                                    , dimnames(original)[[3]]
                                                    )#end expand.grid
                            , MARGIN   = 1
                            , FUN      = paste
                            , collapse = "+"
                            )#end apply
            every   = apply(X = original, FUN = unlist, MARGIN = c(1,4))
            every   = aperm(a = every   , perm = c(2,1,3))
            dimnames(every)[[2]] = concat
         }else{
            every   = eft[[var.vname]]$ts[,,y.sel,,drop=FALSE]
            every   = apply(X = every, FUN = mean, MARGIN = c(1,2,4),na.rm=TRUE)
         }#end if
         #---------------------------------------------------------------------------------#


         #------ Get dimensions. ----------------------------------------------------------#
         r.dim   = dim(every)[-1]      ; names(r.dim  ) = c("real","season")
         r.dname = dimnames(every)[-1] ; names(r.dname) = c("real","season")
         r.seq   = mapply( FUN = sequence, nvec = r.dim, SIMPLIFY = FALSE)
         o.dim   = dim(every)[-1]      ; names(o.dim  ) = c("real","season")
         o.dname = dimnames(every)[-1] ; names(o.dname) = c("real","season")
         o.seq   = mapply( FUN = sequence, nvec = o.dim, SIMPLIFY = FALSE)
         #---------------------------------------------------------------------------------#



         #----- Find the panel information. -----------------------------------------------#
         p.type         = simul$dim.type == "panel"
         p.name         = "panel"
         p.ind          = simul$panel$index  [g.sel]
         p.key          = simul$key[g.sel,p.type]
         p.level        = list(simul$panel$level)
         p.dim          = simul$panel$n.level
         p.seq          = sequence(p.dim)
         names(p.level) = "panel"
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Retrieve the scenarios and the default for each scenario.                  #
         #---------------------------------------------------------------------------------#
         s.type         = simul$dim.type == "scenario"
         s.name         = names(simul$dim[s.type])
         s.dim          = simul$dim[s.type]
         s.dnames       = simul$panel$level
         s.value        = simul$value[g.sel,s.type]
         s.key          = simul$key[g.sel,s.type]
         s.level        = split(x=simul$scenario$level,f=simul$scenario$idxcol)
         names(s.level) = names(simul$dim)[s.type]
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the array names.                                                      #
         #---------------------------------------------------------------------------------#
         s.default = split(simul$default[s.type],index(simul$default[s.type]))
         s.ind     = data.frame( apply(X=simul$scenario$index,MARGIN=2,FUN="[",g.sel)
                               , stringsAsFactors = FALSE
                               )#end data.frame
         a.dim     = sapply(X=c(p.dim,s.dim  ,o.dim  ),FUN=c)
         a.dnames  = c(p.level,s.level,o.dname)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Build the lists to split the array by simulation and panel.                #
         #---------------------------------------------------------------------------------#
         if (! is.list(p.ind)) p.ind = list(p.ind)
         if (! is.list(s.ind)) s.ind = list(s.ind)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Split the stuff into panels, and create the scenario maps.                 #
         #---------------------------------------------------------------------------------#
         this       = array(data=every,dim=a.dim,dimnames=a.dnames)
         n.a.dim    = length(a.dim)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Find the dimensions of each scenario.                                        #
         #---------------------------------------------------------------------------------#
         sc1   = scenario[[2]]
         n.sc1 = length(sc1$key)
         sc2   = scenario[[1]]
         n.sc2 = length(sc2$key)
         n.lax = n.sc2 + 2
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Set up the box plot colour.                                                 #
         #---------------------------------------------------------------------------------#
         bp.colour = rep( x     = c("transparent",sc2$colour,"transparent")
                        , times = n.sc1
                        )#end rep
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the x based on the 1st scenario (n.sc1) and the 2nd scenario          #
         # (n.lax).                                                                        #
         #---------------------------------------------------------------------------------#
         xlimit = c(0,n.sc1*n.lax) + 0.5
         xat    = seq(from=1 + 0.5*(n.lax-1), to=n.sc1*n.lax, by=n.lax)
         xgrid  = seq(from=0.5, to=(n.sc1+1)*n.lax - 0.5*(n.lax-1),by=n.lax)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the y limits.                                                         #
         #---------------------------------------------------------------------------------#
         ylimit.season = pretty.xylim( u       = every[,,-n.season]
                                     , fracexp = 0.0
                                     , is.log  = var.plog
                                     )#end pretty.xylim
         ylimit.year   = pretty.xylim( u       = every[,, n.season]
                                     , fracexp = 0.0
                                     , is.log  = var.plog
                                     )#end pretty.xylim
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Plot the data by season.                                                   #
         #---------------------------------------------------------------------------------#
         for (e in 1:n.season){
            #------------------------------------------------------------------------------#
            #     Keep the current season only.                                            #
            #------------------------------------------------------------------------------#
            bye      = n.a.dim
            inds     = as.matrix(expand.grid(mapply(FUN=sequence,a.dim[-bye])))
            inds     = cbind(inds,e)
            this.var = array(this[inds],dim=a.dim[-bye],dimnames=a.dnames[-bye])
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Set the annotation and some auxiliary variables.                        #
            #------------------------------------------------------------------------------#
            letitre = paste("Scenario comparision - ",var.desc,"\n",season.desc[e]
                           ," - ",global.desc,sep="")
            #------------------------------------------------------------------------------#



            #----- Axis labels. -----------------------------------------------------------#
            lex     = sc1$alabel
            ley     = desc.unit(desc=var.desc,unit=var.unit)
            #------------------------------------------------------------------------------#


            #----- Y limit. ---------------------------------------------------------------#
            ylimit = if(e == n.season){ylimit.year}else{ylimit.season}
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all output formats.                                            #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Get the path and create file name. ----------------------------------#
               if ( e == n.season){
                  now.outpath  = outpath[[o]]$global[[g]]$panelbox.year
                  fichier      = paste(now.outpath,"/panelbox-year-"
                                      ,var.vname,"-",global.suffix,".",outform[o],sep="")
               }else{
                  now.outpath  = outpath[[o]]$global[[g]]$panelbox.season[e]
                  fichier      = paste(now.outpath,"/panelbox-",season.suffix[e],"-"
                                      ,var.vname,"-",global.suffix,".",outform[o],sep="")
               }#end if
               #---------------------------------------------------------------------------#


               #----- Open the file. ------------------------------------------------------#
               if (outform[o] == "x11"){
                  X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=wsize$width*depth
                     ,height=wsize$height*depth,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=wsize$width,height=wsize$height
                            ,pointsize=ptsz,paper=wsize$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=wsize$width
                     ,height=wsize$height,pointsize=ptsz,paper=wsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Save the margins to avoid losing the data. --------------------------#
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Split the plotting window.                                           #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma = c(0.2,3,4.5,0))
               layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
                     , heights = c(rep(5/lo.panel$nrow,lo.panel$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Plot the legend.                                                      #
               #---------------------------------------------------------------------------#
               par(mar=c(0.2,4.4,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "bottom"
                     , inset  = 0.0
                     , legend = if (any(sc2$parse)){
                                   parse(text=sc2$legend)
                                }else{
                                   sc2$legend
                                }#end if
                     , fill   = sc2$colour
                     , border = foreground
                     , ncol   = min(n.sc2,max(3,pretty.box(n=n.sc2)$ncol))
                     , title  = eval(substitute(expression(bold(a)),list(a=sc2$alabel)))
                     , cex    = 0.9
                     , xpd    = TRUE
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over all panels, and plot the data for this season.             #
               #---------------------------------------------------------------------------#
               loop.panel = sequence(max(1,simul$panel$n.level   ))
               for (p in loop.panel){

                  #------------------------------------------------------------------------#
                  #     Find out where is this box going, and set up axes and margins.     #
                  #------------------------------------------------------------------------#
                  left    = (p %% lo.panel$ncol) == 1
                  right   = (p %% lo.panel$ncol) == 0
                  top     = p <= lo.panel$ncol
                  bottom  = p > (lo.panel$nrow - 1) * lo.panel$ncol
                  mar.now = c(1 + 4*bottom,1 + 3*left,1 + 4*top,1 + 3*right) + 0.1
                  #------------------------------------------------------------------------#



                  #----- Make the sub-title. ----------------------------------------------#
                  lesub=simul$panel$title[p]
                  #------------------------------------------------------------------------#



                  #----- Create the temporary box plot list. ------------------------------#
                  pp.this = this.var[p,,,]
                  d1 = which(mapply(identical,dimnames(pp.this),MoreArgs=list(sc1$key)))
                  d2 = which(mapply(identical,dimnames(pp.this),MoreArgs=list(sc2$key)))
                  empty   = array( data     = NA
                                 , dim      = dim(pp.this)[-d2]
                                 , dimnames = dimnames(pp.this)[-d2]
                                 )#end array
                  pp.this = abind(aaa=empty,pp.this,zzz=empty,along=d2)
                  ai      = arrayInd(sequence(length(pp.this)),.dim=dim(pp.this))
                  bp.plot = split(x=pp.this,f=list(ai[,d2],ai[,d1]))
                  #------------------------------------------------------------------------#


                  #----- Set the window. --------------------------------------------------#
                  par(mar = mar.now)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  if (bottom) axis(side=1,at=xat,labels=sc1$label)
                  if (left)   axis(side=2)
                  box()
                  title(main=lesub,cex.main=cex.ptsz,line=0.5)
                  if (plotgrid) abline(h=axTicks(2),col=grid.colour,lty="solid")
                  if (plotgrid) abline(v=xgrid,col=grey.fg,lty="solid")

                  #----- Add the box plot, without the x axis. ----------------------------#
                  boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                         ,yaxt="n")
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               gtitle( main      = letitre
                     , xlab      = lex
                     , ylab      = ley
                     , off.xlab  = 1/6
                     , line.xlab = 4.1
                     , line.ylab = 2.6
                     , cex.main  = 1.1*cex.ptsz
                     )#end gtitle
               # par(las=0)
               # mtext(side=1,text=lex,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
               # mtext(side=2,text=ley,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
               # mtext(side=3,text=letitre,outer=TRUE,cex=1.1*cex.ptsz,font=2)
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               clean.tmp()
               #---------------------------------------------------------------------------#
            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end for (e in 1:n.season)
         #---------------------------------------------------------------------------------#
      }#end if
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
   }#end for (v in 1:npanel.box)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#





   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #    Check whether to create XYZ panel comparisons for sites.                           #
   #---------------------------------------------------------------------------------------#
   if (n.scenario == 2){
      if (is.na(plot.panelxyz)){
         loop.x = sequence(npanel.xvar)[1]
         loop.y = sequence(npanel.yvar)[1]
         loop.z = sequence(npanel.zvar)[1]
      }else if (plot.panelxyz){
         loop.x = sequence(npanel.xvar)
         loop.y = sequence(npanel.yvar)
         loop.z = sequence(npanel.zvar)
      }else{
         loop.x = integer(0)
         loop.y = integer(0)
         loop.z = integer(0)
      }#end if
   }else{
      loop.x = integer(0)
      loop.y = integer(0)
      loop.z = integer(0)
   }#end if
   cat  ("         ~ Plotting the parameter space comparison...","\n")




   #---------------------------------------------------------------------------------------#
   #      Retrieve the scenarios and the default for each scenario.                        #
   #---------------------------------------------------------------------------------------#
   s.type         = simul$dim.type == "scenario"
   s.name         = names(simul$dim[s.type])
   s.dim          = simul$dim[s.type]
   s.dnames       = simul$panel$level
   s.value        = simul$value[g.sel,s.type]
   s.key          = simul$key[g.sel,s.type]
   s.level        = split(x=simul$scenario$level,f=simul$scenario$idxcol)
   names(s.level) = names(simul$dim)[s.type]
   s.default      = split(simul$default[s.type],index(simul$default[s.type]))
   s.ind          = data.frame( apply(X=simul$scenario$index,MARGIN=2,FUN="[",g.sel)
                              , stringsAsFactors = FALSE
                              )#end data.frame
   #---------------------------------------------------------------------------------------#



   #----- Find the panel information. -----------------------------------------------------#
   p.type         = simul$dim.type == "panel"
   p.name         = "panel"
   p.ind          = simul$panel$index  [g.sel]
   p.key          = simul$key[g.sel,p.type]
   p.level        = list(simul$panel$level)
   p.dim          = simul$panel$n.level
   p.seq          = sequence(p.dim)
   names(p.dim)   = "panel"
   names(p.level) = "panel"
   p.title        = gsub(pattern="Phenology: ",replacement="",x=simul$panel$title)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over scenarios.                                                              #
   #---------------------------------------------------------------------------------------#
   for (s in loop.allscen){
      if (s == 0){
         idx.s         = 1
         odx.s         = match(s,loop.allscen)
         out.desc      = simul$global$title  [g]
         out.suffix    = simul$global$level  [g]
         scen.headline = "All scenarios"
      }else{
         idx.s         = match(s,loop.allscen)
         odx.s         = idx.s
         out.desc      = paste(simul$scenario$title[s],simul$global$title  [g],sep=" - ")
         out.suffix    = paste(simul$scenario$level[s],simul$global$level  [g],sep="-")
         scen.headline = simul$scenario$title[s]
      }#end if

      #------------------------------------------------------------------------------------#
      #    Find the column with the plot information.                                      #
      #------------------------------------------------------------------------------------#
      if (simul$scenario$n.level == 0){
         now = scenario[[1]]
      }else if (s == 0){
         now = scenario[-simul$scenario$idxcol[1]][[1]]
      }else{
         now = scenario[-simul$scenario$idxcol[s]][[1]]
      }#end if
      n.now = length(now$key)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Select the runs that belong to this panel.                                     #
      #------------------------------------------------------------------------------------#
      if (simul$scenario$n.level == 0 | s == 0){
         s.sel = rep(TRUE,times=n.gsel)
      }else{
         s.sel = simul$scenario$idxoff[g.sel,simul$scenario$idxcol[s]] == s
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over Y variables.                                                         #
      #------------------------------------------------------------------------------------#
      for (y in loop.y){
         y.vname    = panel.xyz$yvar$vname   [y]
         y.desc     = panel.xyz$yvar$desc    [y]
         y.unit     = panel.xyz$yvar$unit    [y]
         y.plog     = panel.xyz$yvar$plog    [y]
         y.plog.dbh = panel.xyz$yvar$plog.dbh[y]
         y.leg      = panel.xyz$yvar$leg     [y]

         cat  ("           = Y variable: ",y.desc," scenario ",scen.headline,"...","\n")

         #----- Grab the time series for Y. -----------------------------------------------#
         y.ts    = eft[[y.vname]]$ts      [s.sel,,y.sel,n.season]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Grab the dimensions.                                                       #
         #---------------------------------------------------------------------------------#
         r.dim   = dim(y.ts)[-1]      ; names(r.dim  ) = c("real","year")
         r.dname = dimnames(y.ts)[-1] ; names(r.dname) = c("real","year")
         r.seq   = mapply( FUN = sequence, nvec = r.dim, SIMPLIFY = FALSE)
         o.dim   = dim(y.ts)[-1]      ; names(o.dim  ) = c("real","year")
         o.dname = dimnames(y.ts)[-1] ; names(o.dname) = c("real","year")
         o.seq   = mapply( FUN = sequence, nvec = o.dim, SIMPLIFY = FALSE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Find the array names.  For the special case in which we keep all           #
         # scenarios, we must re-arrange arrays differently.                               #
         #---------------------------------------------------------------------------------#
         if (s == 0){
            a.dim     = sapply(X=c(p.dim,s.dim[idx.s],s.dim[-idx.s],o.dim  ),FUN=c)
            a.dnames  = c(p.level,s.level[idx.s],s.level[-idx.s],o.dname)
            n.a.dim   = length(a.dim)
            d.1st     = 3
         }else{
            a.dim     = sapply(X=c(p.dim,s.dim[-idx.s],o.dim  ),FUN=c)
            a.dnames  = c(p.level,s.level[-idx.s],o.dname)
            n.a.dim   = length(a.dim)
            d.1st     = 2
         }#end if
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #      Split the stuff into panels, and create the scenario maps.                 #
         #---------------------------------------------------------------------------------#
         y.ts     = aperm( a    = array(data=y.ts,dim=a.dim,dimnames=a.dnames)
                         , perm = c(d.1st,sequence(n.a.dim)[-d.1st])
                         )#end aperm
         #---------------------------------------------------------------------------------#


         #----- Create indices to split the arrays into lists. ----------------------------#
         ai.ts    = arrayInd(sequence(length(y.ts   )),.dim=dim(y.ts   ))
         #---------------------------------------------------------------------------------#



         #----- Find the format of the points for the plots. ------------------------------#
         pch.ts    = array( data = now$pch[arrayInd(1:length(y.ts   )
                                                   ,.dim=dim(y.ts   ))[,1]]
                          , dim  = dim(y.ts)
                          )#end array
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Loop over X variables.                                                      #
         #---------------------------------------------------------------------------------#
         for (x in loop.x){
            x.vname    = panel.xyz$xvar$vname   [x]
            x.desc     = panel.xyz$xvar$desc    [x]
            x.unit     = panel.xyz$xvar$unit    [x]
            x.plog     = panel.xyz$xvar$plog    [x]
            x.plog.dbh = panel.xyz$xvar$plog.dbh[x]
            x.leg      = panel.xyz$xvar$leg     [x]
            x.pft      = panel.xyz$xvar$pftvar  [x]
            x.dbh      = panel.xyz$xvar$dbhvar  [x]


            #------------------------------------------------------------------------------#
            #      Split the stuff into panels, and create the scenario maps.              #
            #------------------------------------------------------------------------------#
            x.ts     = eft[[x.vname]]$ts     [s.sel,,y.sel,n.season]
            x.ts     = aperm( a    = array(data=x.ts,dim=a.dim,dimnames=a.dnames)
                            , perm = c(d.1st,sequence(n.a.dim)[-d.1st])
                            )#end aperm
            #---------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Plot only when x and y are different.                                    #
            #------------------------------------------------------------------------------#
            if (x.vname != y.vname){
               #----- Build the log info. -------------------------------------------------#
               xy.plog=paste(ifelse(x.plog,"x",""),ifelse(y.plog,"y",""),sep="")
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over Z variables.                                                #
               #---------------------------------------------------------------------------#
               for (z in loop.z){
                  z.vname    = panel.xyz$zvar$vname         [z]
                  z.desc     = panel.xyz$zvar$desc          [z]
                  z.unit     = panel.xyz$zvar$unit          [z]
                  z.plog     = panel.xyz$zvar$plog          [z]
                  z.plog.dbh = panel.xyz$zvar$plog.dbh      [z]
                  z.plog.xyz = panel.xyz$zvar$plog.xyz      [z]
                  z.pft      = panel.xyz$zvar$pftvar        [z]
                  z.dbh      = panel.xyz$zvar$dbhvar        [z]
                  z.cscheme  = get(panel.xyz$zvar$col.scheme[z])


                  #------------------------------------------------------------------------#
                  #     Plot only when x and y are different.                              #
                  #------------------------------------------------------------------------#
                  if ( ! z.vname %in% c(x.vname,y.vname)){
                     #---------------------------------------------------------------------#
                     #      Split the stuff into panels, and create the scenario maps.     #
                     #---------------------------------------------------------------------#
                     z.ts       = eft[[z.vname]]$ts   [s.sel,,y.sel,n.season]
                     z.ts     = aperm( a    = array(data=z.ts,dim=a.dim,dimnames=a.dnames)
                                     , perm = c(d.1st,sequence(n.a.dim)[-d.1st])
                                     )#end aperm
                     #---------------------------------------------------------------------#
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                     #     Split the arrays into lists.                                    #
                     #---------------------------------------------------------------------#
                     x.list    = split(x=x.ts  ,f=ai.ts[,2])
                     y.list    = split(x=y.ts  ,f=ai.ts[,2])
                     z.list    = split(x=z.ts  ,f=ai.ts[,2])
                     pch.list  = split(x=pch.ts,f=ai.ts[,2])
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #      Fix the titles.                                                #
                     #---------------------------------------------------------------------#
                     letitre   = paste(z.desc,"\n",out.desc,sep="")
                     lex       = desc.unit(desc=x.desc,unit=x.unit)
                     ley       = desc.unit(desc=y.desc,unit=y.unit)
                     lacle     = desc.unit(desc=NULL  ,unit=z.unit)
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Plot the parameter space by season.                             #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Get the path. ----------------------------------------------#
                        now.outpath  = outpath[[o]]$global[[g]]$panel.xyz[odx.s,y]
                        #------------------------------------------------------------------#



                        #----- Open file or display. --------------------------------------#
                        fichier = paste(now.outpath,"/y-",y.vname,"_x-",x.vname
                                       ,"_z-",z.vname,"-",out.suffix,"-panelcomp."
                                       ,outform[o],sep="")
                        if (outform[o] == "x11"){
                           X11(width=size$width,height=size$height,pointsize=ptsz)
                        }else if(outform[o] == "png"){
                           png(filename=fichier,width=size$width*depth
                              ,height=size$height*depth,pointsize=ptsz,res=depth)
                        }else if(outform[o] == "eps"){
                           postscript(file=fichier,width=size$width,height=size$height
                                     ,pointsize=ptsz,paper=size$paper)
                        }else if(outform[o] == "pdf"){
                           pdf(file=fichier,onefile=FALSE,width=size$width
                              ,height=size$height,pointsize=ptsz,paper=size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Plot the boxes.                                             #
                        #------------------------------------------------------------------#
                        leg.ncol    = min(3,pretty.box(n.now)$ncol)
                        leg.letitre = expression(bold("Simulation"))
                        leg.legend  = if (any(now$parse)){
                                         parse(text=now$legend)
                                      }else{
                                         now$legend
                                      }#end if

                        par(par.user)
                        xyz.plot( x              = x.list
                                , y              = y.list
                                , z              = z.list
                                , fixed.xlim     = TRUE
                                , fixed.ylim     = TRUE
                                , xy.log         = xy.plog
                                , pch            = pch.list
                                , cex            = 1.2
                                , nlevels        = n.colourbar
                                , colour.palette = z.cscheme
                                , xyz.sub        = p.title
                                , xyz.title      = list( main      = letitre
                                                       , xlab      = lex
                                                       , ylab      = ley
                                                       , line.xlab = 4.1
                                                       , line.ylab = 2.6
                                                       , cex.main  = 1.1 * cex.ptsz
                                                       )#end list
                                , xyz.more       = list(grid=list(col=grid.colour
                                                                 ,lty="solid"))
                                , key.log        = z.plog.xyz
                                , key.title      = list(main=lacle,cex.main=0.8*cex.main)
                                , xyz.legend     = list( x      = "bottom"
                                                       , inset  = 0.0
                                                       , legend = leg.legend
                                                       , pch    = now$pch
                                                       , col    = grey.fg
                                                       , border = foreground
                                                       , bg     = background
                                                       , ncol   = leg.ncol
                                                       , title  = leg.letitre
                                                       , cex    = 1.0
                                                       , xpd    = TRUE
                                                       )#end legend
                                )#end xyz.plot
                        #------------------------------------------------------------------#



                        #----- Close the device. ------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        clean.tmp()
                        #------------------------------------------------------------------#
                     }#end for
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                  }#end if (! z.vname %in% c(x.vname,y.vname))
                  #------------------------------------------------------------------------#
               }#end for (z in loop.z)
               #---------------------------------------------------------------------------#
            }#end if (x.vname != y.vname)
            #------------------------------------------------------------------------------#
         }#end for (x in loop.x)
         #---------------------------------------------------------------------------------#
      }#end for (y in loop.y)
      #------------------------------------------------------------------------------------#
   }#end for (s in loop.allscen)
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
   #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

   rm(eft)
}#end for (g in loop.global)
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Save status file so we know it is complete.                                          #
#------------------------------------------------------------------------------------------#
cat (" + Saving the simulation status to ",basename(rdata.status),"...","\n")
cat (" + Saving the simulation status to ",basename(rdata.status),"...","\n")
status = data.frame( rj = n.total, total = n.total, complete = rj == n.total)
write.table( x         = status
           , file      = rdata.status
           , append    = FALSE
           , quote     = FALSE
           , row.names = FALSE
           , col.names = TRUE
           )#end write.table
#==========================================================================================#
#==========================================================================================#
