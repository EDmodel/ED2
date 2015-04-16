#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
here          = getwd()                     # Current directory
srcdir        = "/prj/prjidfca/marcosl/Util/Rsc" # Script directory
stext.default = "stext16"                   # Default soil texture
drain.default = "r-160"                     # Default rainfall
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
                     , colour  = c("#FED164","#FDA531","#FF5308","#D90B00","#6E0500")
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
plot.vec.axes  = FALSE                 # Plot vector axes? (TRUE|FALSE)
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


#----- Make sure that the all directories exist. ------------------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
pcadbhroot = rep(x=NA_character_,times=nout)
pcaallroot = rep(x=NA_character_,times=nout)
for (o in sequence(nout)){
   onow.main     = file.path(outroot,outform[o])
   pcadbhroot[o] = file.path(onow.main,"pcadbh_year")
   pcaallroot[o] = file.path(onow.main,"pcaall_year")
   if (! outform[o] %in% c("quartz","x11")){
      if (! file.exists(onow.main    )) dir.create(onow.main    )
      if (! file.exists(pcadbhroot[o])) dir.create(pcadbhroot[o])
      if (! file.exists(pcaallroot[o])) dir.create(pcaallroot[o])
   }#end if (! outform[o] %in% c("quartz","x11"))
   #---------------------------------------------------------------------------------------#
}#end for (o in sequence(nout))
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
if (! file.exists(rdata.path)) dir.create(rdata.path)
rdata.siminfo = file.path( rdata.path
                         , paste( comp.prefix,stext.default,drain.default,"SimInfo.RData"
                                , sep = "_"
                                )#end paste
                         )#end file.path
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#     Here we read or reload the simulations.                                              #
#------------------------------------------------------------------------------------------#
if (file.exists(rdata.siminfo)){
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
      loop.rdims = seq_along(rdims)
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
#      Time to plot.  Here we must loop over all scenarios, nested within the panels,      #
# which are nested within the global.                                                      #
#------------------------------------------------------------------------------------------#
cat(" + Processing data for:","\n",sep="")
if (is.logical(use.global) && all(use.global)){
   loop.global = sequence(max(1,simul$global$n.level  ))
}else{
   loop.global = use.global
}#end if
loop.panel     = sequence(max(1,simul$panel$n.level))
loop.scenario  = which(simul$scenario$level %in% simul$default)
loop.allscen   = c(loop.scenario,0)


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




   #----- Grab data from previously loaded variable. --------------------------------------#
   cat("   - Retrieving data from file :",basename(rdata.global),"...","\n")
   load(rdata.global)
   eft = get(simul$global$level[g])
   rm(list=simul$global$level[g])
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make some arrays with the right dimensions.                                       #
   #---------------------------------------------------------------------------------------#
   n.real         = simul$n.real
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
   #---------------------------------------------------------------------------------------#




   #----- Create the list of variables. ---------------------------------------------------#
   addpftdbh = c("agb","nmon.wdef","rain","smpot","water.deficit")
   for (a in seq_along(addpftdbh)){
      v.now                 = addpftdbh[a]
      eft[[v.now]]$tspft    = tspft.array
      eft[[v.now]]$tspftdbh = tspftdbh.array
      for (f in sequence(n.pft)){
         eft[[v.now]]$tspft[,,,,f] = eft[[v.now]]$ts
         for (d in sequence(n.dbh)){
            eft[[v.now]]$tspftdbh[,,,,f,d] = eft[[v.now]]$ts
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over all panels.                                                             #
   #---------------------------------------------------------------------------------------#
   for (p in loop.panel){
      #------------------------------------------------------------------------------------#
      #     Select the runs that belong to this panel.                                     #
      #------------------------------------------------------------------------------------#
      if (simul$panel$n.level == 0){
         p.sel = rep(TRUE,times=n.gsel)
         cat("   - Single panel","\n",sep="")
      }else{
         p.sel = simul$panel$index[g.sel] == p
         cat("   - Panel: ",p," - ",simul$panel$title[p],"\n",sep="")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over scenarios.                                                           #
      #------------------------------------------------------------------------------------#
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
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Find the column with the plot information.                                   #
         #---------------------------------------------------------------------------------#
         if (simul$scenario$n.level == 0){
            now = scenario[[1]]
         }else if (s == 0){
            now = scenario[-simul$scenario$idxcol[1]][[1]]
         }else{
            now = scenario[-simul$scenario$idxcol[s]][[1]]
         }#end if
         n.now = length(now$key)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Select the runs that belong to this scenario.                               #
         #---------------------------------------------------------------------------------#
         if (simul$scenario$n.level == 0 | s == 0){
            s.sel = rep(TRUE,times=n.gsel)
         }else{
            s.sel = simul$scenario$idxoff[g.sel,simul$scenario$idxcol[s]] == s
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Selection flag combines both this panel and this scenario.                  #
         #---------------------------------------------------------------------------------#
         sel = p.sel & s.sel
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Build the total description, and the total suffix.                         #
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
            #---- Select the current scenario based on the index. -------------------------#
            if (s == 0){
               scen.title.now = "All scenarios"
               scen.level.now = "allscen"
            }else{
               scen.title.now = simul$scenario$title[s]
               scen.level.now = simul$scenario$level[s]
            }#end if
            #------------------------------------------------------------------------------#

            if (simul$panel$n.level == 0){
               if (simul$global$n.level == 0){
                  #----- 2-D scenarios but no panels or global variables. -----------------#
                  out.desc   = scen.title.now
                  out.suffix = scen.level.now
                  #------------------------------------------------------------------------#
               }else{
                  #----- 2-D scenarios and global variables. ------------------------------#
                  out.desc   = paste(scen.title.now
                                    ,simul$global$title  [g]
                                    ,sep=" - ")
                  out.suffix = paste(scen.level.now
                                    ,simul$global$level  [g]
                                    ,sep="-")
                  #------------------------------------------------------------------------#
               }#end if (simul$global$n.level == 0)
               #---------------------------------------------------------------------------#
            }else{
               if (simul$global$n.level == 0){
                  #----- 2-D scenarios and panel variables. -------------------------------#
                  out.desc   = paste(scen.title.now
                                    ,simul$panel$title   [p]
                                    ,sep=" - ")
                  out.suffix = paste(scen.level.now
                                    ,simul$panel$level   [p]
                                    ,sep="-")
                  #------------------------------------------------------------------------#
               }else{
                  #----- 2-D scenarios, panel, and global variables. ----------------------#
                  out.desc   = paste(paste(scen.title.now
                                          ,simul$panel$title   [p]
                                          ,sep=" - ")
                                    ,simul$global$title  [g]
                                    ,sep="\n")
                  out.suffix = paste(scen.level.now
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
         #     Loop over all DBH classes.                                                  #
         #---------------------------------------------------------------------------------#
         pca.dbh = list()
         for (d in sequence(n.dbh+1)){
            if ( d == n.dbh+1 ){
               cat("     * Community-wide...","\n")
            }else{
               cat("     * DBH class: ",dbh.desc[d],"...","\n")
            }#end if

            #------------------------------------------------------------------------------#
            #     Loop over all explanatory variables.                                     #
            #------------------------------------------------------------------------------#
            for (ev in sequence(npca.explain)){
               exp.vname = pca.explain$vname [ev]
               exp.desc  = pca.explain$desc  [ev]
               exp.unit  = pca.explain$unit  [ev]

               #----- Load the data. ------------------------------------------------------#
               if ( d == n.dbh+1 ){
                  now      = c(eft[[exp.vname]]$ts      [sel,,y.sel,n.season])
               }else{
                  now      = c(eft[[exp.vname]]$tspftdbh[sel,,y.sel,n.season,d,n.pft])
               }#end if
               now[!is.finite(now)] = NA
               mean.now = mean(now,na.rm=TRUE)
               sdev.now = sd  (now,na.rm=TRUE)
               if (! is.finite(mean.now)) mean.now = NA
               if (! is.finite(sdev.now) | sdev.now == 0.) sdev.now = NA
               now = ( now - mean.now ) / sdev.now
               now[! is.finite(now)] = NA
               if ( ev == 1 ){
                  datum              = data.frame( x = now )
                  names(datum)       = exp.vname
               }else{
                  datum[[exp.vname]] = now
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (ev in sequence(npca.explain))
            keep  = rowSums(is.na(datum)) == 0
            datum = datum[keep,]
            #------------------------------------------------------------------------------#



            #----- Calculate the principal component. -------------------------------------#
            pca.dbh[[d]] = prcomp(datum)
            #------------------------------------------------------------------------------#
         }#end for (d in sequence(n.dbh+1))
         #---------------------------------------------------------------------------------#











         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #     Plot all PCAs with labels near the arrows.                                  #
         #---------------------------------------------------------------------------------#
         cat("     * Plot PCA by DBH...","\n")
         for (o in sequence(nout)){
            #----- Open file or display. --------------------------------------------------#
            fichier = file.path( pcadbhroot[o]
                               , paste( "noleg-pcadbh-year-",out.suffix
                                      , ".",outform[o],sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=ssize$width,height=ssize$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=ssize$width,height=ssize$height,pointsize=ptsz
                         ,paper=ssize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                  ,pointsize=ptsz,paper=ssize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split the panel. -------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2,0))
            layout(mat=lo.dbh$mat)
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #     Loop over DBH classes, and plot the PCA analysis.                        #
            #------------------------------------------------------------------------------#
            for (d in sequence(n.dbh)){
               #----- Load PCA. -----------------------------------------------------------#
               pca.now  = pca.dbh[[d]]
               summ.pca = summary(pca.now)
               scores   = pca.now$x
               rota     = pca.now$rotation
               npca     = nrow(scores)
               lam      = pca.now$sdev * sqrt(npca)
               pca.pts  = t(t(scores) / lam)
               pca.vec  = t(t(rota  ) * lam)
               explain  = round(100*summ.pca$importance[2,],1)
               sig.vec  = sqrt(rota[,1]^2+rota[,2]^2) >= 0.1
               lwd.vec  = ifelse(sig.vec,      2,       2)
               lty.vec  = ifelse(sig.vec,"solid","dashed")
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #   Make axes.                                                              #
               #---------------------------------------------------------------------------#
               pvar.explained = paste("Component",sequence(npca)
                                     ," - ", sprintf("%.1f",explain),"%",sep="")
               #---------------------------------------------------------------------------#



               #------ Find limits. -------------------------------------------------------#
               pxlimit  = max(abs(pca.pts[,1])) * c(-1,1)
               pylimit  = max(abs(pca.pts[,2])) * c(-1,1)
               vxlimit  = max(abs(pca.vec[,1])) * c(-2.0,2.0)
               vylimit  = max(abs(pca.vec[,2])) * c(-2.0,2.0)
               #---------------------------------------------------------------------------#


               #----- Plot the PCA points. ------------------------------------------------#
               par(mar=c(4.1,4.1,5.1,2.1))
               plot.new()
               plot.window(xlim=pxlimit,ylim=pylimit)
               abline(h=0,v=0,col=foreground,lty="solid",lwd=2)
               axis(side=1)
               axis(side=2)
               title( main = paste("DBH Class: ",dbh.desc[d],sep="")
                    , xlab = pvar.explained[1]
                    , ylab = pvar.explained[2]
                    )#end title
               points(x=pca.pts[,1],y=pca.pts[,2],col=washed.mg,cex=0.5,pch=16)
               #---------------------------------------------------------------------------#


               #----- Plot the PCA vectors. -----------------------------------------------#
               plot.window(xlim=vxlimit,ylim=vylimit)
               if (plot.vec.axes){
                  axis(side=3,col.ticks=firebrick.fg,col.axis=firebrick.mg)
                  axis(side=4,col.ticks=firebrick.fg,col.axis=firebrick.mg)
               }#end if
               box()
               for (u in sequence(npca.explain)){
                  arrows(x0=0,y0=0,x1=pca.vec[u,1],y1=pca.vec[u,2]
                        ,col=pca.explain$colour[u],length=0.10
                        ,lwd=lwd.vec[u],lty=lty.vec[u])
               }#end for
               for (u in sequence(npca.explain)){
                  rot = 180*atan2(pca.vec[u,2],pca.vec[u,1])/pi
                  if (rot >  90 & rot <=  180) rot = rot + 180
                  if (rot < -90 & rot >= -180) rot = rot + 180
                  mult = 1.50 + runif(n=1,min=-0.20,max=0.20)
                  text  (x=mult*pca.vec[u,1],y=mult*pca.vec[u,2]
                        ,labels=parse(text=pca.explain$short[u]),col=pca.explain$colour[u]
                        ,cex=0.9,srt=rot)
               }#end for
               #---------------------------------------------------------------------------#
            }#end for (d in sequence(n.dbh))
            #------------------------------------------------------------------------------#



            #----- Main title. ------------------------------------------------------------#
            gtitle( main = out.desc )
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for for (o in sequence(nout))
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#











         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #     Plot all PCAs with labels near the arrows.                                  #
         #---------------------------------------------------------------------------------#
         cat("     * Plot community-wide PCA...","\n")
         for (o in sequence(nout)){
            #----- Open file or display. --------------------------------------------------#
            fichier = file.path( pcaallroot[o]
                               , paste( "noleg-pcaall-year-",out.suffix
                                      , ".",outform[o],sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=ssize$width,height=ssize$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=ssize$width,height=ssize$height,pointsize=ptsz
                         ,paper=ssize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                  ,pointsize=ptsz,paper=ssize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Load PCA. --------------------------------------------------------------#
            pca.now  = pca.dbh[[n.dbh+1]]
            summ.pca = summary(pca.now)
            scores   = pca.now$x
            rota     = pca.now$rotation
            npca     = nrow(scores)
            lam      = pca.now$sdev * sqrt(npca)
            pca.pts  = t(t(scores) / lam)
            pca.vec  = t(t(rota  ) * lam)
            explain  = round(100*summ.pca$importance[2,],1)
            sig.vec  = sqrt(rota[,1]^2+rota[,2]^2) >= 0.1
            lwd.vec  = ifelse(sig.vec,      2,       2)
            lty.vec  = ifelse(sig.vec,"solid","dashed")
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #   Make axes.                                                                 #
            #------------------------------------------------------------------------------#
            pvar.explained = paste("Component",sequence(npca)
                                  ," - ", sprintf("%.1f",explain),"%",sep="")
            #------------------------------------------------------------------------------#



            #------ Find limits. ----------------------------------------------------------#
            pxlimit  = max(abs(pca.pts[,1])) * c(-1,1)
            pylimit  = max(abs(pca.pts[,2])) * c(-1,1)
            vxlimit  = max(abs(pca.vec[,1])) * c(-2.0,2.0)
            vylimit  = max(abs(pca.vec[,2])) * c(-2.0,2.0)
            #------------------------------------------------------------------------------#


            #----- Plot the PCA points. ---------------------------------------------------#
            par(par.user)
            plot.new()
            plot.window(xlim=pxlimit,ylim=pylimit)
            abline(h=0,v=0,col=foreground,lty="solid",lwd=2)
            axis(side=1)
            axis(side=2)
            title( main = out.desc
                 , xlab = pvar.explained[1]
                 , ylab = pvar.explained[2]
                 )#end title
            points(x=pca.pts[,1],y=pca.pts[,2],col=washed.mg,cex=0.5,pch=16)
            #------------------------------------------------------------------------------#





            #----- Plot the PCA vectors. --------------------------------------------------#
            plot.window(xlim=vxlimit,ylim=vylimit)
            if (plot.vec.axes){
               axis(side=3,col.ticks=firebrick.fg,col.axis=firebrick.mg)
               axis(side=4,col.ticks=firebrick.fg,col.axis=firebrick.mg)
            }#end if
            box()
            for (u in sequence(npca.explain)){
               arrows(x0=0,y0=0,x1=pca.vec[u,1],y1=pca.vec[u,2]
                     ,col=pca.explain$colour[u],length=0.10
                     ,lwd=lwd.vec[u],lty=lty.vec[u])
            }#end for
            for (u in sequence(npca.explain)){
               rot = 180*atan2(pca.vec[u,2],pca.vec[u,1])/pi
               if (rot >  90 & rot <=  180) rot = rot + 180
               if (rot < -90 & rot >= -180) rot = rot + 180
               mult = 1.50 + runif(n=1,min=-0.20,max=0.20)
               text  (x=mult*pca.vec[u,1],y=mult*pca.vec[u,2]
                     ,labels=parse(text=pca.explain$short[u]),col=pca.explain$colour[u]
                     ,cex=0.9,srt=rot)
            }#end for
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for for (o in sequence(nout))
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#










         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #     Plot all PCAs with legends at the bottom.                                   #
         #---------------------------------------------------------------------------------#
         cat("     * Plot PCA by DBH with legend...","\n")
         for (o in sequence(nout)){
            #----- Open file or display. --------------------------------------------------#
            fichier = file.path( pcadbhroot[o]
                               , paste( "legend-pcadbh-year-",out.suffix
                                      , ".",outform[o],sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=psize$width,height=psize$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height,pointsize=ptsz
                         ,paper=psize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
                  ,pointsize=ptsz,paper=psize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split the panel. -------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2,0))
            layout( mat     = rbind(lo.dbh$mat.off,rep(1,times=lo.dbh$ncol))
                  , heights = c(rep(5/lo.dbh$nrow,times=lo.dbh$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     First, the legend.                                                       #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "center"
                  , inset  = 0.0
                  , legend = pca.explain$lname
                  , fill   = pca.explain$colour
                  , border = foreground
                  , title  = expression(bold("Variables - Annual Means"))
                  , xpd    = TRUE
                  , cex    = 0.9 * cex.ptsz
                  , ncol   = 6 # min(5,pretty.box(npca.explain)$ncol)
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Loop over DBH classes, and plot the PCA analysis.                        #
            #------------------------------------------------------------------------------#
            for (d in sequence(n.dbh)){
               #----- Load PCA. -----------------------------------------------------------#
               pca.now  = pca.dbh[[d]]
               summ.pca = summary(pca.now)
               scores   = pca.now$x
               rota     = pca.now$rotation
               npca     = nrow(scores)
               lam      = pca.now$sdev * sqrt(npca)
               pca.pts  = t(t(scores) / lam)
               pca.vec  = t(t(rota  ) * lam)
               explain  = round(100*summ.pca$importance[2,],1)
               sig.vec  = sqrt(rota[,1]^2+rota[,2]^2) >= 0.1
               lwd.vec  = ifelse(sig.vec,      2,       2)
               lty.vec  = ifelse(sig.vec,"solid","dashed")
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #   Make axes.                                                              #
               #---------------------------------------------------------------------------#
               pvar.explained = paste("Component",sequence(npca)
                                     ," - ", sprintf("%.1f",explain),"%",sep="")
               #---------------------------------------------------------------------------#



               #------ Find limits. -------------------------------------------------------#
               pxlimit  = max(abs(pca.pts[,1])) * c(-1,1)
               pylimit  = max(abs(pca.pts[,2])) * c(-1,1)
               vxlimit  = max(abs(pca.vec[,1])) * c(-1,1)
               vylimit  = max(abs(pca.vec[,2])) * c(-1,1)
               #---------------------------------------------------------------------------#


               #----- Plot the PCA points. ------------------------------------------------#
               par(mar=c(4.1,4.1,5.1,2.1))
               plot.new()
               plot.window(xlim=pxlimit,ylim=pylimit)
               abline(h=0,v=0,col=foreground,lty="solid",lwd=2)
               axis(side=1)
               axis(side=2)
               title( main = paste("DBH Class: ",dbh.desc[d],sep="")
                    , xlab = pvar.explained[1]
                    , ylab = pvar.explained[2]
                    )#end title
               points(x=pca.pts[,1],y=pca.pts[,2],col=washed.bg,cex=0.5,pch=16)
               #---------------------------------------------------------------------------#


               #----- Plot the PCA vectors. -----------------------------------------------#
               plot.window(xlim=vxlimit,ylim=vylimit)
               if (plot.vec.axes){
                  axis(side=3,col.ticks=firebrick.fg,col.axis=firebrick.mg)
                  axis(side=4,col.ticks=firebrick.fg,col.axis=firebrick.mg)
               }#end if
               box()
               for (u in sequence(npca.explain)){
                  arrows(x0=0,y0=0,x1=pca.vec[u,1],y1=pca.vec[u,2]
                        ,col=pca.explain$colour[u],length=0.10
                        ,lwd=lwd.vec[u],lty=lty.vec[u])
               }#end for
               #---------------------------------------------------------------------------#
            }#end for (d in sequence(n.dbh))
            #------------------------------------------------------------------------------#



            #----- Main title. ------------------------------------------------------------#
            gtitle( main = out.desc )
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for for (o in sequence(nout))
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#










         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #     Plot all PCAs with legends at the bottom.                                   #
         #---------------------------------------------------------------------------------#
         cat("     * Plot community-wide PCA with legend...","\n")
         for (o in sequence(nout)){
            #----- Open file or display. --------------------------------------------------#
            fichier = file.path( pcaallroot[o]
                               , paste( "legend-pcaall-year-",out.suffix
                                      , ".",outform[o],sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=psize$width,height=psize$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height,pointsize=ptsz
                         ,paper=psize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
                  ,pointsize=ptsz,paper=psize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split the panel. -------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2,0))
            layout(mat= rbind(2,1),heights=c(5,1))
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     First, the legend.                                                       #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "center"
                  , inset  = 0.0
                  , legend = pca.explain$lname
                  , fill   = pca.explain$colour
                  , border = foreground
                  , title  = expression(bold("Variables - Annual Means"))
                  , xpd    = TRUE
                  , cex    = 0.8 * cex.ptsz
                  , ncol   = 6 # min(5,pretty.box(npca.explain)$ncol)
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Load and plot the PCA analysis.                                          #
            #------------------------------------------------------------------------------#
            pca.now  = pca.dbh[[n.dbh+1]]
            summ.pca = summary(pca.now)
            scores   = pca.now$x
            rota     = pca.now$rotation
            npca     = nrow(scores)
            lam      = pca.now$sdev * sqrt(npca)
            pca.pts  = t(t(scores) / lam)
            pca.vec  = t(t(rota  ) * lam)
            explain  = round(100*summ.pca$importance[2,],1)
            sig.vec  = sqrt(rota[,1]^2+rota[,2]^2) >= 0.1
            lwd.vec  = ifelse(sig.vec,      2,       2)
            lty.vec  = ifelse(sig.vec,"solid","dashed")
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #   Make axes.                                                                 #
            #------------------------------------------------------------------------------#
            pvar.explained = paste("Component",sequence(npca)
                                  ," - ", sprintf("%.1f",explain),"%",sep="")
            #------------------------------------------------------------------------------#



            #------ Find limits. ----------------------------------------------------------#
            pxlimit  = max(abs(pca.pts[,1])) * c(-1,1)
            pylimit  = max(abs(pca.pts[,2])) * c(-1,1)
            vxlimit  = max(abs(pca.vec[,1])) * c(-1,1)
            vylimit  = max(abs(pca.vec[,2])) * c(-1,1)
            #------------------------------------------------------------------------------#


            #----- Plot the PCA points. ---------------------------------------------------#
            par(mar=c(4.1,4.1,5.1,2.1))
            plot.new()
            plot.window(xlim=pxlimit,ylim=pylimit)
            abline(h=0,v=0,col=foreground,lty="solid",lwd=2)
            axis(side=1)
            axis(side=2)
            title( main = out.desc
                 , xlab = pvar.explained[1]
                 , ylab = pvar.explained[2]
                 )#end title
            points(x=pca.pts[,1],y=pca.pts[,2],col=washed.bg,cex=0.5,pch=16)
            #------------------------------------------------------------------------------#


            #----- Plot the PCA vectors. --------------------------------------------------#
            plot.window(xlim=vxlimit,ylim=vylimit)
            if (plot.vec.axes){
               axis(side=3,col.ticks=firebrick.fg,col.axis=firebrick.mg)
               axis(side=4,col.ticks=firebrick.fg,col.axis=firebrick.mg)
            }#end if
            box()
            for (u in sequence(npca.explain)){
               arrows(x0=0,y0=0,x1=pca.vec[u,1],y1=pca.vec[u,2]
                     ,col=pca.explain$colour[u],length=0.10
                     ,lwd=lwd.vec[u],lty=lty.vec[u])
            }#end for
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for for (o in sequence(nout))
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      }#end for (s in loop.allscen)
      #------------------------------------------------------------------------------------#
   }#end for (p in loop.panel)
   #---------------------------------------------------------------------------------------#
}#end for (g in loop.global)
#==========================================================================================#
#==========================================================================================#
