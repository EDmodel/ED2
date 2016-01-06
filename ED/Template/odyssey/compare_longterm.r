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
here    = getwd()                               #   Current directory
srcdir  = "/n/home00/mlongo/util/Rsc"           #   Script directory
ibackground    = 0                              # Make figures compatible to background
                                                # 0 -- white
                                                # 1 -- black
                                                # 2 -- dark grey
#----- Output directory -------------------------------------------------------------------#
outroot = file.path(here,paste0("longterm_comp_ibg",sprintf("%2.2i",ibackground)))
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.hour  = c(FALSE,TRUE)[2]
reload.range = c(FALSE,TRUE)[2]
rdata.path   = file.path(here,"RData_longterm")
rdata.suffix = "longterm_ed22.RData"
#------------------------------------------------------------------------------------------#



#----- Default settings. ------------------------------------------------------------------#
emean.yeara = 2002  # First year
emean.yearz = 2005  # Last year
#------------------------------------------------------------------------------------------#



#----- Default settings. ------------------------------------------------------------------#
pdens.range = c(0.05,1.0)      # First year
pdens.log   = c(FALSE,TRUE)[1] # Plot density function as logarithm?
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# Type of DBH class                                                                        #
#    1 -- Every 10 cm until 100cm; > 100cm                                                 #
#    2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)                                      #
#    3 -- 0-10; 10-35; 35-55; > 55 (cm)                                                    #
#------------------------------------------------------------------------------------------#
idbh.type = 3
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
eort  = "t"
n     = 0
sites = list()
n          = n + 1
sites[[n]] = list( iata = "gyf"
                 , desc = "Paracou"
                 , pch  = 2
                 , col  = "#3B24B3"
                 , fg   = "#160959"
                 , drya = "08/31"
                 , dryz = "11/05"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "s67"
                 , desc = "Santarem km 67"
                 , pch  =  5
                 , col  = "#A3CC52"
                 , fg   = "#4B6614"
                 , drya = "07/13"
                 , dryz = "11/21"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "pdg"
                 , desc = "Pe-de-Gigante"
                 , pch  = 13
                 , col  = "#990F0F"
                 , fg   = "#4D0404"
                 , drya = "04/17"
                 , dryz = "09/26"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "rja"
                 , desc = "Rebio Jaru"
                 , pch  =  1
                 , col  = "#306614"
                 , fg   = "#143305"
                 , drya = "05/13"
                 , dryz = "09/09"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "m34"
                 , desc = "Manaus K34"
                 , pch  =  6
                 , col  = "#2996CC"
                 , fg   = "#0A4766"
                 , drya = "07/04"
                 , dryz = "10/01"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "pnz"
                 , desc = "Petrolina"
                 , pch  =  4
                 , col  = "#B49ED2"
                 , fg   = "#7D6E93"
                 , drya = "03/18"
                 , dryz = "01/16"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "ban"
                 , desc = "Bananal"
                 , pch  =  8
                 , col  = "#F5C858"
                 , fg   = "#AB8C3D"
                 , drya = "05/10"
                 , dryz = "09/29"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "bsb"
                 , desc = "Brasilia"
                 , pch  =  9
                 , col  = "#E65C17"
                 , fg   = "#732A06"
                 , drya = "04/19"
                 , dryz = "10/02"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "nat"
                 , desc = "Natal"
                 , pch  =  0
                 , col  = "#00F3FB"
                 , fg   = "#00AAAF"
                 , drya = "08/19"
                 , dryz = "02/07"
                 )#end list
n          = n + 1
sites[[n]] = list( iata = "tnf"
                 , desc = "Tapajos National Forest"
                 , pch  =  5
                 , col  = "#A3CC52"
                 , fg   = "#4B6614"
                 , drya = "07/13"
                 , dryz = "11/21"
                 )#end list
use.sites  = "tnf"
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.struct  = list( name        = c("ihrz00_irad01_ccislp010"
                                   ,"ihrz02_irad01_ccislp010"
                                   ,"ihrz02_irad01_ccislp030"
                                   )#end c
                  , desc        = c("Horizontal OFF"
                                   ,"Original CCI correction"
                                   ,"Enhanced CCI correction"
                                   )#end c
                  , verbose     = c("Horizontal OFF"
                                   ,"Original CCI correction"
                                   ,"Enhanced CCI correction"
                                   )#end c
                  , colour     = c("#3B24B3","#2996CC","#E65C17")
                  , fgcol      = c("#160959","#0A4766","#732A06")
                  , age.interp = c(     TRUE,     TRUE,     TRUE)
                  )#end list
#----- List the default simulation. -------------------------------------------------------#
sim.default = 1
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

depth             = 300                  # PNG/TIFF resolution, in pixels per inch
paper             = "square"             # Paper size, to define the plot shape
ptsz              = 17                   # Font size.
n.dens            = 512                  # Number of density points. 
dens.min          = 0.00001              # Minimum density (relative to maximum)
pft.use           = c(2,3,4,18)          # PFTs to use
ncolours.xyz      = 96                   # # colours for xyz plots
f.leg             = 1/6                  # Factor to expand plot devices
#----- Maximum absolute value for Bowen ratio. --------------------------------------------#
bmn = -0.09
bmx =  1.99
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Switch controls to plot only the needed ones.                                       #
#------------------------------------------------------------------------------------------#
plot.ts.emean   = c(FALSE,TRUE)[2]
plot.ts.ymean   = c(FALSE,TRUE)[2]
plot.ts.pft     = c(FALSE,TRUE)[2]
plot.ts.dbh     = c(FALSE,TRUE)[2]
plot.mm.pft     = c(FALSE,TRUE)[2]
plot.mm.dbh     = c(FALSE,TRUE)[2]
plot.pdf.patch  = c(FALSE,TRUE)[2]
plot.xyz.patch  = c(FALSE,TRUE)[2]
plot.zm.patch   = c(FALSE,TRUE)[2]
plot.ym.patch   = c(FALSE,TRUE)[2]
plot.ym.theme   = c(FALSE,TRUE)[2]
col.dryseason   = "papayawhip"
col.ust.altern  = "firebrick4"
col.ust.default = "deeppink"
slz.cscheme     = "visible"
patch.aggr      = c("age","lorey")[1]
#------------------------------------------------------------------------------------------#



#----- Treefall disturbance rate. ---------------------------------------------------------#
treefall.default = 0.0111
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
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Eddy flux comparisons.                                                               #
#------------------------------------------------------------------------------------------#
n = 0
compvar       = list()
n             = n + 1
compvar[[ n]] = list( vnam         = "agb"
                    , desc         = "Above-ground biomass"
                    , unit         = untab$kgcom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = "nplant"
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "ba"
                    , desc         = "Basal area"
                    , unit         = untab$m2om2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = "nplant"
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "lai"
                    , desc         = "Leaf area index"
                    , unit         = untab$m2lom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "gpp"
                    , desc         = "Gross primary productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "npp"
                    , desc         = "Net primary productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "plant.resp"
                    , desc         = "Plant respiration"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "cba"
                    , desc         = "Carbon balance"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "reco"
                    , desc         = "Ecosystem respiration"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "nep"
                    , desc         = "Net Ecosystem Productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "hflxca"
                    , desc         = "Sensible heat flux"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "wflxca"
                    , desc         = "Water vapour flux"
                    , unit         = untab$kgwom2oday
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "transp"
                    , desc         = "Transpiration"
                    , unit         = untab$kgwom2oday
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "bowen"
                    , desc         = "Bowen ratio"
                    , unit         = untab$empty
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "tratio"
                    , desc         = "Transpiration ratio"
                    , unit         = untab$empty
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "ustar"
                    , desc         = "Friction velocity"
                    , unit         = untab$mos
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "rshortup"
                    , desc         = "Upward SW radiation"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "albedo"
                    , desc         = "Albedo"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "rlongup"
                    , desc         = "Upward LW radiation"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "parup"
                    , desc         = "Upward PAR"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "par.gnd"
                    , desc         = "Ground absorption - PAR"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "rshort.gnd"
                    , desc         = "Ground absorption - SW"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.gpp"
                    , desc         = "Leaf GPP"
                    , unit         = untab$kgcom2loyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.gsw"
                    , desc         = "Stomatal conductance"
                    , unit         = untab$kgwom2loday
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.temp"
                    , desc         = "Leaf temperature"
                    , unit         = untab$degC
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.vpd"
                    , desc         = "Leaf VPD"
                    , unit         = untab$hpa
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "par.leaf"
                    , desc         = "Leaf Absorption - PAR"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "par.leaf.beam"
                    , desc         = "Leaf Absorption - Direct PAR"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "par.leaf.diff"
                    , desc         = "Leaf Absorption - Diffuse PAR"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.par"
                    , desc         = "Norm. Leaf Absorption - PAR"
                    , unit         = untab$umolom2los
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.par.beam"
                    , desc         = "Norm. Leaf Absorption - Direct PAR"
                    , unit         = untab$umolom2los
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "leaf.par.diff"
                    , desc         = "Norm. Leaf Absorption - Diffuse PAR"
                    , unit         = untab$umolom2los
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "assim.light"
                    , desc         = "Light-limited Assimilation"
                    , unit         = untab$umolom2los
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "assim.rubp"
                    , desc         = "RuBP-limited Assimilation"
                    , unit         = untab$umolom2los
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "assim.co2"
                    , desc         = "CO2-limited Assimilation"
                    , unit         = untab$umolom2los
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "fast.soil.c"
                    , desc         = "Fast soil carbon"
                    , unit         = untab$kgcom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "struct.soil.c"
                    , desc         = "Structural soil carbon"
                    , unit         = untab$kgcom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "slow.soil.c"
                    , desc         = "Slow soil carbon"
                    , unit         = untab$kgcom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "tot.soil.c"
                    , desc         = "Soil carbon"
                    , unit         = untab$kgcom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "can.depth"
                    , desc         = "Mean canopy height"
                    , unit         = untab$m
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "can.area"
                    , desc         = "Mean canopy area"
                    , unit         = untab$empty
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "wood.dens"
                    , desc         = "Mean wood density"
                    , unit         = untab$gocm3
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "firemort"
                    , desc         = "Fire mortality"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "ncbmort"
                    , desc         = "Density-dependent mortality"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "dimort"
                    , desc         = "Density independent mortality"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "mort"
                    , desc         = "Mortality rate"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Set the themes.                                                                      #
#------------------------------------------------------------------------------------------#
n             = 0
ym.theme      = list()
n             = n + 1
ym.theme[[n]] = list( vnam   = c("fast.soil.c","slow.soil.c","struct.soil.c","tot.soil.c")
                    , desc   = c("Fast"       ,"Slow"       ,"Structural"   ,"Total"     )
                    , colour = c("#A3CC52"    ,"#3B24B3"    ,"#E65C17"      ,"#6B6B6B"   )
                    , lwd    = c(2.5          ,2.5          ,2.5            ,2.5         )
                    , ylog   = FALSE
                    , prefix = "soil_carbon"
                    , title  = "Soil Carbon"
                    , unit   = untab$kgcom2
                    )#end list
n             = n + 1
ym.theme[[n]] = list( vnam   = c("can.depth"   )
                    , desc   = c("Canopy depth")
                    , colour = c("#A3CC52"     )
                    , lwd    = c(2.5           )
                    , ylog   = FALSE
                    , prefix = "can_depth"
                    , title  = "Canopy height"
                    , unit   = untab$m
                    )#end list
n             = n + 1
ym.theme[[n]] = list( vnam   = c("can.area"    )
                    , desc   = c("Canopy area" )
                    , colour = c("#A3CC52"     )
                    , lwd    = c(2.5           )
                    , ylog   = FALSE
                    , prefix = "can_area"
                    , title  = "Canopy area"
                    , unit   = untab$empty
                    )#end list
n             = n + 1
ym.theme[[n]] = list( vnam   = c("wood.dens"   )
                    , desc   = c("Wood density")
                    , colour = c("#E65C17"     )
                    , lwd    = c(2.5           )
                    , ylog   = FALSE
                    , prefix = "wood_dens"
                    , title  = "Wood density"
                    , unit   = untab$gocm3
                    )#end list
n             = n + 1
ym.theme[[n]] = list( vnam   = c("firemort","ncbmort"  ,"dimort"       ,"mort"      )
                    , desc   = c("Fire"
                                ,"Neg. carbon balance"
                                ,"Density independent"
                                ,"Total"
                                )#end c
                    , colour = c("#E65C17" ,"#3B24B3"  ,"#A3CC52"      ,"#6B6B6B"   )
                    , lwd    = c(2.5        ,2.5            ,2.5         )
                    , ylog   = FALSE
                    , prefix = "mortality"
                    , title  = "Mortality rates"
                    , unit   = untab$pcpopoyr
                    )#end list
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#

#----- Create directory with RData. -------------------------------------------------------#
if (! file.exists(rdata.path)) dir.create(rdata.path)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Combine all structures into a consistent list.                                       #
#------------------------------------------------------------------------------------------#
n.sim       = length(sim.struct$name)
#----- Simulation keys. -------------------------------------------------------------------#
simul.key   = sim.struct$name
#----- Description. -----------------------------------------------------------------------#
simleg.key  = sim.struct$desc
#---- Create the colours and line type for legend. ----------------------------------------#
simcol.key  = sim.struct$colour
simfgc.key  = sim.struct$fgcol
simlty.key  = rep("solid",times=n.sim)
simcex.key  = rep(2.0    ,times=n.sim)
simain.key  = sim.struct$age.interp
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Dump the information to a list.                                                      #
#------------------------------------------------------------------------------------------#
simul       = data.frame( name             = simul.key
                        , desc             = simleg.key
                        , colour           = simcol.key
                        , fgcol            = simfgc.key
                        , lty              = simlty.key
                        , cex              = simcex.key
                        , age.interp       = simain.key
                        , stringsAsFactors = FALSE
                        )#end data.frame
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Convert sites to data frame, and keep only those to be run.                         #
#------------------------------------------------------------------------------------------#
sites        = list.2.data.frame(sites)
if (is.logical(use.sites)){
   use.sites = which(rep(use.sites,times=nrow(sites))[sequence(nrow(sites))])
}else if (is.character(use.sites)){
   use.sites = match(use.sites,sites$iata)
   use.sites = use.sites[! is.na(use.sites)]
}#end if
if (length(use.sites) == 0){
   stop ("You must specify at least one valid site")
}else{
   sites = sites[use.sites,]
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      List the keys for all dimensions.                                                   #
#------------------------------------------------------------------------------------------#
compvar.key  = list.2.data.frame(compvar)$vnam
compvar.pch  = list.2.data.frame(compvar)$pch
season.key   = season.list
moment.key   = c("mean","variance","skewness","kurtosis")
moment.desc  = c("Mean","Variance","Skewness","Kurtosis")
#------------------------------------------------------------------------------------------#




#----- Set the various dimensions associated with variables, simulations, and sites. ------#
nsites     = length(sites$iata )
nsimul     = length(simul.key  )
ncompvar   = length(compvar.key)
nym.theme  = length(ym.theme   )
nseason    = length(season.key )
nmoment    = length(moment.key )
#------------------------------------------------------------------------------------------#



#----- Suffix for seasons (so it stays in order). -----------------------------------------#
season.suffix = paste(sprintf("%2.2i",sequence(nseason)),tolower(season.key),sep="-")
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     List for size plots.                                                                 #
#------------------------------------------------------------------------------------------#
dbh.use     = sequence(ndbh)
season.use  = sequence(nseason-1)
npft.use    = length(pft.use)
ndbh.use    = length(dbh.use)
nseason.use = length(season.use)
#------------------------------------------------------------------------------------------#



#----- Find the best set up for plotting all seasons in the same plot. --------------------#
lo.box    = pretty.box(n=nseason-1)
lo.simul  = pretty.box(n=nsimul,byrow=FALSE)
lo.pft    = pretty.box(n=npft.use)
lo.dbh    = pretty.box(n=ndbh.use)
lo.site   = pretty.box(n=nsites)
lo.season = pretty.box(n=nseason.use)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (1. - f.leg)
size0   = plotsize(proje=FALSE,paper=paper)
height0 = size0$height
width0  = size0$width
xsize   = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
#----- Scale for monthly and annual means. ------------------------------------------------#
eywidth  = width0  * n.sim
eyheight = height0
eysize   = plotsize( proje     = FALSE
                   , stdheight = eyheight
                   , stdwidth  = eywidth
                   , extendfc  = "lat"
                   , extfactor = f.ext
                   )#end plotsize
#----- Scale for multi-simulation plots. --------------------------------------------------#
swidth  = width0  * max(1.,lo.simul$ncol / lo.simul$nrow)
sheight = height0 * max(1.,lo.simul$nrow / lo.simul$ncol)
ssize   = plotsize( proje     = FALSE
                  , stdheight = sheight
                  , stdwidth  = swidth
                  , extendfc  = "lat"
                  , extfactor = f.ext
                  )#end plotsize
#----- Scale for multi-PFT plots. ---------------------------------------------------------#
fwidth  = width0  * max(1.,lo.pft$ncol / lo.pft$nrow)
fheight = height0 * max(1.,lo.pft$nrow / lo.pft$ncol)
fsize   = plotsize( proje     = FALSE
                  , stdheight = fheight
                  , stdwidth  = fwidth
                  , extendfc  = "lat"
                  , extfactor = f.ext
                  )#end plotsize
#----- Scale for multi-DBH plots. ---------------------------------------------------------#
dwidth  = width0  * max(1.,lo.dbh$ncol / lo.dbh$nrow)
dheight = height0 * max(1.,lo.dbh$nrow / lo.dbh$ncol)
dsize   = plotsize( proje     = FALSE
                  , stdheight = dheight
                  , stdwidth  = dwidth
                  , extendfc  = "lat"
                  , extfactor = f.ext
                  )#end plotsize
#----- Scale for multi-DBH plots. ---------------------------------------------------------#
zwidth  = width0  * max(1.,lo.season$ncol / lo.season$nrow)
zheight = height0 * max(1.,lo.season$nrow / lo.season$ncol)
zsize   = plotsize( proje     = FALSE
                  , stdheight = zheight
                  , stdwidth  = zwidth
                  , extendfc  = "lat"
                  , extfactor = f.ext
                  )#end plotsize
#----- Scale for multi-site plots. --------------------------------------------------------#
pwidth  = width0  * max(1.,lo.site$ncol  / lo.site$nrow )
pheight = height0 * max(1.,lo.site$nrow  / lo.site$ncol )
psize   = plotsize( proje     = FALSE
                  , stdheight = sheight
                  , stdwidth  = swidth
                  , extendfc  = "lat"
                  , extfactor = f.ext
                  )#end plotsize
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
out = list()
for (o in sequence(nout)){
   is.figure   = ! outform[o] %in% c("quartz","x11")
   this.form   = outform[o]


   #---- Main path for this output format. ------------------------------------------------#
   o.form = list()
   o.form$main = file.path(outroot,this.form)
   if (is.figure && ! file.exists(o.form$main)) dir.create(o.form$main)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of monthly means.                                #
   #---------------------------------------------------------------------------------------#
   if (plot.ts.emean){
      o.ts.emean = file.path(o.form$main,"ts_emean")
      if (is.figure){
         if (! file.exists(o.ts.emean)) dir.create(o.ts.emean)
      }#end if
      o.form$ts.emean = o.ts.emean
   }#end if (plot.ts.emean)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of monthly means.                                #
   #---------------------------------------------------------------------------------------#
   if (plot.ts.ymean){
      o.ts.ymean = file.path(o.form$main,"ts_ymean")
      if (is.figure){
         if (! file.exists(o.ts.ymean)) dir.create(o.ts.ymean)
      }#end if
      o.form$ts.ymean = o.ts.ymean
   }#end if (plot.ts.emean)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of annual means.                                 #
   #---------------------------------------------------------------------------------------#
   if (plot.ts.pft){
      o.main   = file.path(o.form$main,"ts_pft")
      o.ts.pft = list( main      = o.main
                     , default   = file.path(o.main,"default"  )
                     , variables = file.path(o.main,"variables")
                     )#end list
      if (is.figure){
         if (! file.exists(o.ts.pft$main     )) dir.create(o.ts.pft$main     )
         if (! file.exists(o.ts.pft$default  )) dir.create(o.ts.pft$default  )
         if (! file.exists(o.ts.pft$variables)) dir.create(o.ts.pft$variables)
      }#end if
      o.form$ts.pft = o.ts.pft
   }#end if (plot.ts.pft)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of annual means.                                 #
   #---------------------------------------------------------------------------------------#
   if (plot.ts.dbh){
      o.main   = file.path(o.form$main,"ts_dbh")
      o.ts.dbh = list( main      = o.main
                     , default   = file.path(o.main,"default"  )
                     , variables = file.path(o.main,"variables")
                     )#end list
      if (is.figure){
         if (! file.exists(o.ts.dbh$main     )) dir.create(o.ts.dbh$main     )
         if (! file.exists(o.ts.dbh$default  )) dir.create(o.ts.dbh$default  )
         if (! file.exists(o.ts.dbh$variables)) dir.create(o.ts.dbh$variables)
      }#end if
      o.form$ts.dbh = o.ts.dbh
   }#end if (plot.ts.dbh)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the monthly means by PFT.                                        #
   #---------------------------------------------------------------------------------------#
   if (plot.mm.pft){
      o.main   = file.path(o.form$main,"mm_pft")
      o.mm.pft = list( main      = o.main
                     , default   = file.path(o.main,"default"  )
                     , variables = file.path(o.main,"variables")
                     )#end list
      if (is.figure){
         if (! file.exists(o.mm.pft$main     )) dir.create(o.mm.pft$main     )
         if (! file.exists(o.mm.pft$default  )) dir.create(o.mm.pft$default  )
         if (! file.exists(o.mm.pft$variables)) dir.create(o.mm.pft$variables)
      }#end if
      o.form$mm.pft = o.mm.pft
   }#end if (plot.mm.pft)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the monthly means by DBH.                                        #
   #---------------------------------------------------------------------------------------#
   if (plot.mm.dbh){
      o.main   = file.path(o.form$main,"mm_dbh")
      o.mm.dbh = list( main      = o.main
                     , default   = file.path(o.main,"default"  )
                     , variables = file.path(o.main,"variables")
                     )#end list
      if (is.figure){
         if (! file.exists(o.mm.dbh$main     )) dir.create(o.mm.dbh$main     )
         if (! file.exists(o.mm.dbh$default  )) dir.create(o.mm.dbh$default  )
         if (! file.exists(o.mm.dbh$variables)) dir.create(o.mm.dbh$variables)
      }#end if
      o.form$mm.dbh = o.mm.dbh
   }#end if (plot.mm.dbh)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the ym by time of year and patch.                               #
   #---------------------------------------------------------------------------------------#
   if (plot.ym.patch){
      o.main  = file.path(o.form$main,"ym_patch")
      o.ym.patch = list( main      = o.main
                       , default   = file.path(o.main,"default")
                       , variables = file.path(o.main,"variables")
                       )#end list
      if (is.figure){
         if (! file.exists(o.ym.patch$main     )) dir.create(o.ym.patch$main     )
         if (! file.exists(o.ym.patch$default  )) dir.create(o.ym.patch$default  )
         if (! file.exists(o.ym.patch$variables)) dir.create(o.ym.patch$variables)
      }#end if
      o.form$ym.patch = o.ym.patch
   }#end if (plot.ym.patch)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the zm by season and patch property.                             #
   #---------------------------------------------------------------------------------------#
   if (plot.zm.patch){
      o.main  = file.path(o.form$main,"zm_patch")
      o.zm.patch = list( main      = o.main
                       , default   = file.path(o.main,"default")
                       , variables = file.path(o.main,"variables")
                       )#end list
      if (is.figure){
         if (! file.exists(o.zm.patch$main     )) dir.create(o.zm.patch$main     )
         if (! file.exists(o.zm.patch$default  )) dir.create(o.zm.patch$default  )
         if (! file.exists(o.zm.patch$variables)) dir.create(o.zm.patch$variables)
      }#end if
      o.form$zm.patch = o.zm.patch
   }#end if (plot.zm.patch)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the PDF by time of year and patch.                               #
   #---------------------------------------------------------------------------------------#
   if (plot.pdf.patch){
      o.main  = file.path(o.form$main,"pdf_patch")
      o.pdf.patch = list( main      = o.main
                        , default   = file.path(o.main,"default"  )
                        , variables = file.path(o.main,"variables")
                        )#end list
      if (is.figure){
         if (! file.exists(o.pdf.patch$main     )) dir.create(o.pdf.patch$main     )
         if (! file.exists(o.pdf.patch$default  )) dir.create(o.pdf.patch$default  )
         if (! file.exists(o.pdf.patch$variables)) dir.create(o.pdf.patch$variables)
      }#end if
      o.form$pdf.patch = o.pdf.patch
   }#end if (plot.pdf.patch)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the mean annual cycle by age.                                    #
   #---------------------------------------------------------------------------------------#
   if (plot.xyz.patch){
      o.main  = file.path(o.form$main,"xyz_patch")
      o.xyz.patch = list( main      = o.main
                    , default   = file.path(o.main,"default")
                    , variables = file.path(o.main,"variables")
                    )#end list
      if (is.figure){
         if (! file.exists(o.xyz.patch$main     )) dir.create(o.xyz.patch$main     )
         if (! file.exists(o.xyz.patch$default  )) dir.create(o.xyz.patch$default  )
         if (! file.exists(o.xyz.patch$variables)) dir.create(o.xyz.patch$variables)
      }#end if
      o.form$xyz.patch = o.xyz.patch
   }#end if (plot.xyz.patch)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the theme plots.                                                 #
   #---------------------------------------------------------------------------------------#
   if (plot.ym.theme){
      o.main     = file.path(o.form$main,"ym_theme")
      o.ym.theme = list( main      = o.main
                       , default   = file.path(o.main,"default")
                       , variables = file.path(o.main,"variables")
                       )#end list
      if (is.figure){
         if (! file.exists(o.ym.theme$main     )) dir.create(o.ym.theme$main     )
         if (! file.exists(o.ym.theme$default  )) dir.create(o.ym.theme$default  )
         if (! file.exists(o.ym.theme$variables)) dir.create(o.ym.theme$variables)
      }#end if
      o.form$ym.theme = o.ym.theme
   }#end if (plot.ym.theme)
   #---------------------------------------------------------------------------------------#


   #----- Save the full list to the main path list. ---------------------------------------#
   out[[this.form]] = o.form
   #---------------------------------------------------------------------------------------#
}#end for (o in 1:nout)
#------------------------------------------------------------------------------------------#







#------------------------------------------------------------------------------------------#
#      Retrieve all data that already exists.                                              #
#------------------------------------------------------------------------------------------#
res        = list()
loop.sites = integer(length=0)
for (p in sequence(nsites)){
   iata = sites$iata[p]
   #----- Find file name. -----------------------------------------------------------------#
   rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
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
      loop.sites = c(loop.sites,p)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Loop over all sites and simulations to determine the age range.                      #
#------------------------------------------------------------------------------------------#
#----- Find file name. --------------------------------------------------------------------#
rdata.range = file.path(rdata.path,paste("var_range",rdata.suffix,sep="_"))
if (reload.range && file.exists(rdata.range)){
   #----- Reload ranges. ------------------------------------------------------------------#
   cat0("   - Load the ranges...")
   dummy = load(file=rdata.range)
   #---------------------------------------------------------------------------------------#
}else{
   #----- Find ranges. --------------------------------------------------------------------#
   cat0("   - Find the ranges...")
   age.range   = array(data=NA,dim=c(2,nsites,nsimul))
   lorey.range = array(data=NA,dim=c(2,nsites,nsimul))
   var.range   = array(data=NA,dim=c(2,ncompvar,nsites,nsimul))
   for (p in sequence(nsites)){
      iata          = sites$iata[p]
      longname      = sites$desc[p]
      cat0("    > ",sites$desc[p],"...")
      for (s in sequence(nsimul)){
         cat0("      # Simulation: ",simul$desc[s],"...")

         #----- Load hourly averages. -----------------------------------------------------#
         ans.name = paste0(eort,iata,"_",simul$name[s])
         ans.path = file.path(here,ans.name)
         ans.file = file.path(ans.path,"rdata_month",paste0(ans.name,".RData"))
         load(ans.file)
         #---------------------------------------------------------------------------------#



         #----- Create time stamp for annual and monthly means. ---------------------------#
         nymean  = emean.yearz - emean.yeara + 1
         toyear  = rep(seq(from=emean.yeara,to=emean.yearz,by=1),each=12)
         tomonth = chron( paste( month = rep( sequence(12), times=nymean   )
                               , day   = rep(            1, times=nymean*12)
                               , year  = toyear
                               , sep   = "/"
                               )#end paste
                        )#end chron
         nemean  = tomonth
         #---------------------------------------------------------------------------------#



         #----- Grab all relevant times and . ---------------------------------------------#
         esel       = datum$when %wr% tomonth
         emean.loop = match(datum$when[esel],tomonth)
         agepa      = datum$patch$age
         loreypa    = datum$patch$can.depth
         #---------------------------------------------------------------------------------#



         #----- Update ranges. ------------------------------------------------------------#
         for (e in emean.loop){
            now   = tomonth[e]
            emap  = match(tomonth[e],datum$when)
            mm    = nummonths(now)
            yyyy  = numyears (now)
            stamp = paste0("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm))
            if (simul$age.interp[s]){
               agenow            = pmax(1/6,agepa[[stamp]])
               loreynow          = loreypa[[stamp]]
               age.range  [,p,s] = range(c(age.range  [,p,s],agenow  ),finite=TRUE)
               lorey.range[,p,s] = range(c(lorey.range[,p,s],loreynow),finite=TRUE)
            }#end if


            for (v in sequence(ncompvar)){
               #----- Load information. ---------------------------------------------------#
               this.compvar      = compvar[[v]]
               this.vnam         = this.compvar$vnam
               this.patchvar     = this.compvar$patchvar
               #---------------------------------------------------------------------------#



               #----- Several variables require special handling, check these first. ------#
               if (this.vnam %in% "tot.soil.c" && this.patchvar){
                  #----- tot.soil.c is a combination of variables. ------------------------#
                  vnow = ( datum$patch [["fast.soil.c"  ]][[stamp]]
                         + datum$patch [["struct.soil.c"]][[stamp]]
                         + datum$patch [["slow.soil.c"  ]][[stamp]]
                         )#end vnow
                  #------------------------------------------------------------------------#
               }else if (this.vnam %in% "tot.soil.c"){
                  #----- tot.soil.c is a combination of variables. ------------------------#
                  vnow = with(datum$emean, fast.soil.c  [emap]
                                         + slow.soil.c  [emap]
                                         + struct.soil.c[emap]
                             )#end with
                  #------------------------------------------------------------------------#
               }else if(this.vnam %in% "firemort"){
                  #------ Fire mortality: convert exponential rate to fraction rate. ------#
                  fire    = pmax(0.,datum$lu$dist[emap,3,3]-treefall.default,na.rm=TRUE)
                  vnow    = pmax(0.1,100. * ( 1.0 - exp( - fire)))
                  #------------------------------------------------------------------------#
               }else if(this.vnam %in% "dimort"){
                  #------ Dens. Independent rate: remove fire and make it fraction rate. --#
                  fire    = pmax(0.,datum$lu$dist[emap,3,3]-treefall.default)
                  dimort  = pmax(0.,datum$szpft$dimort[emap,5,18]-fire)
                  vnow    = pmax(0.1,100. * ( 1.0 - exp( - dimort )))
                  #------------------------------------------------------------------------#
               }else if(this.vnam %in% c("mort","ncbmort")){
                  #------ Other mortalities: make them fraction rate. ---------------------#
                  nowmort = pmax(0.0,datum$szpft[[this.vnam]][emap,5,18])
                  vnow    = pmax(0.1,100. * ( 1.0 - exp( - nowmort)))
                  #------------------------------------------------------------------------#
               }else if(this.vnam %in% "bowen"){
                  #------ Bowen ratio: use the patch-level variable. ----------------------#
                  hflxca  = datum$patch$hflxca [[stamp]]
                  qwflxca = datum$patch$qwflxca[[stamp]]
                  vnow    = ifelse(qwflxca %!=% 0, hflxca/qwflxca, NA)
               }else if(this.vnam %in% "tratio"){
                  #------ Bowen ratio: use the patch-level variable. ----------------------#
                  transp  = datum$patch$transp[[stamp]]
                  wflxca  = datum$patch$wflxca[[stamp]]
                  vnow    = ifelse(wflxca %!=% 0, transp/wflxca, NA)
               }else if(this.patchvar){
                  #------ Use patch-level info if available... ----------------------------#
                  vnow = datum$patch [[this.vnam]][[stamp]]
                  #------------------------------------------------------------------------#
               }else{
                  #------ ... otherwise use polygon averaged monthly means... -------------#
                  vnow = datum$emean[[this.vnam]][emap]
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#


               #----- Make sure variables are either finite or NA, and update range. ------#
               vnow = ifelse(is.finite(vnow),vnow,NA)
               if (any(is.finite(vnow))){
                  var.range[,v,p,s] = range(c(var.range[,v,p,s],vnow),finite=TRUE)
               }#end if
               rm(vnow)
               #---------------------------------------------------------------------------#
            }#end for (v in sequence(ncompvar))
            #------------------------------------------------------------------------------#
            rm(now,mm,yyyy,stamp,agenow)
         }#end for (e in emean.loop)
         #---------------------------------------------------------------------------------#
         rm(esel,tomonth,nemean,emean.loop,agepa,loreypa)
      }#end for (s in sequence(nsimul))
      #------------------------------------------------------------------------------------#
   }#end for (p in loop.sites)
   #---------------------------------------------------------------------------------------#


   #----- Make sure the range is either finite or NA. -------------------------------------#
   var.range = ifelse(is.finite(var.range),var.range,NA)
   dimnames(var.range) = list(c("min","max"),compvar.key,sites$iata,simul$name)
   #---------------------------------------------------------------------------------------#


   #----- Save range. ---------------------------------------------------------------------#
   cat0("   - Finding the ranges...")
   dummy = save(list=c("lorey.range","age.range","var.range"),file=rdata.range)
   #---------------------------------------------------------------------------------------#
}#end if (reload.range && file.exists(rdata.range))
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Define age classes.                                                                  #
#------------------------------------------------------------------------------------------#
age.at     = unique(pretty.log(age.range,n=30,forcelog=TRUE))
lnage.at   = log(age.at)
n.age.at   = length(age.at)
lorey.at   = unique(pretty(lorey.range,n=30))
n.lorey.at = length(lorey.at)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Loop over all sites, variables, and simulations, to prepare the data for the        #
# model comparison.  We want to use only times with actual measurements, so we will        #
# discard model results from times with no observation so all derived quantities have      #
# the same number of defined points (so if measurements are biased towards daytime, the    #
# model will also be equally biased).                                                      #
#------------------------------------------------------------------------------------------#
if (length(loop.sites) != 0) cat0(" + Processing missing hourly data...")
for (p in loop.sites){
   #----- Get the basic information. ------------------------------------------------------#
   iata          = sites$iata[p]
   im            = match(iata,poilist$iata)
   this          = list()
   this$short    = poilist$short   [im]
   this$longname = sites$desc      [ p]
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
   for (s in sequence(nsimul)){
      cat0("      # Simulation: ",simul$desc[s],"...")

      #----- Load hourly averages. --------------------------------------------------------#
      ans.name = paste0(eort,iata,"_",simul$name[s])
      ans.path = file.path(here,ans.name)
      ans.file = file.path(ans.path,"rdata_month",paste0(ans.name,".RData"))
      load(ans.file)
      #------------------------------------------------------------------------------------#





      #----- Create some variables to describe season and time of the day. ----------------#
      model = list()



      #----- Create time stamp for annual and monthly means. ------------------------------#
      nymean        = emean.yearz - emean.yeara + 1
      emean.toyear  = rep(seq(from=emean.yeara,to=emean.yearz,by=1),each=12)
      model$tomonth = chron( paste( month = rep( sequence(12), times=nymean   )
                                  , day   = rep(            1, times=nymean*12)
                                  , year  = emean.toyear
                                  , sep   = "/"
                                  )#end paste
                           )#end chron
      model$toyear  = unique(numyears(model$tomonth))
      nemean        = length(model$tomonth)
      #------------------------------------------------------------------------------------#



      #----- Grab all relevant times and . ------------------------------------------------#
      esel       = datum$when %wr% model$tomonth
      emean.loop = match(datum$when[esel],model$tomonth)
      agepa      = datum$patch$age
      loreypa    = datum$patch$can.depth
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #      Load all variables, interpolate them and make the table.                      #
      #------------------------------------------------------------------------------------#
      cat0("       ~ Load variables...")
      for (v in sequence(ncompvar)){
         #----- Load information. ---------------------------------------------------------#
         this.compvar   = compvar[[v]]
         this.vnam      = this.compvar$vnam
         this.desc      = this.compvar$desc
         this.unit      = this.compvar$unit
         this.szpftvar  = this.compvar$szpftvar
         this.patchvar  = this.compvar$patchvar
         this.cohortvar = this.compvar$cohortvar
         this.scalevar  = this.compvar$scalevar
         cat0("         > ",this.desc,"...")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Grab all formats of this variable (it is fine if some of them are null).    #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% "tot.soil.c"){
            emean = with(datum$emean,fast.soil.c+struct.soil.c+slow.soil.c)
            szpft = NULL
            if (this.patchvar){
               patch = with(datum$patch,mapply( FUN      =  function(fast,slow,struct){
                                                               out = fast + slow + struct
                                                               return(out)
                                                            }#end function
                                              , fast     = fast.soil.c
                                              , slow     = slow.soil.c
                                              , struct   = struct.soil.c
                                              , SIMPLIFY = FALSE
                                              )#end mapply
                           )#end with
            }else{
               patch = NULL
            }#end if (this.patchvar)
            #------------------------------------------------------------------------------#
         }else if (this.vnam %in% "ba"){
            emean = datum$szpft[[this.vnam]][,5,18]
            szpft = datum$szpft[[this.vnam]]
            patch = datum$patch[[this.vnam]]
         }else if (this.vnam %in% "firemort"){
            emean = pmax(0.0,datum$lu$dist[,3,3] - treefall.default)
            szpft = array( data     = emean
                         , dim      = dim(datum$szpft$mort)
                         , dimnames = dimnames(datum$szpft$mort)
                         )#end array
            patch = NULL
         }else if (this.vnam %in% "dimort"){
            emean.fire = pmax(0.0,datum$lu$dist[,3,3] - treefall.default)
            szpft.fire = array(data=emean.fire,dim=dim(datum$szpft$mort))
            emean      = pmax(0.0,datum$szpft[[this.vnam]][,5,18]-emean.fire)
            szpft      = ifelse( datum$szpft[[this.vnam]] > szpft.fire
                               , datum$szpft[[this.vnam]] - szpft.fire
                               , 0.0
                               )#end ifelse
            patch      = NULL
         }else if (this.vnam %in% c("mort","ncbmort")){
            emean = datum$szpft[[this.vnam]][,5,18]
            szpft = datum$szpft[[this.vnam]]
            patch = NULL
         }else if (this.vnam %in% "tot.dist"){
            emean = 100. * pmax(datum$lu$dist[,3,3],0.0,na.rm=TRUE)
            szpft = NULL
            patch = NULL
         }else if (this.vnam %in% "bowen"){
            emean = with(datum$emean,ifelse(qwflxca %!=% 0,hflxca/qwflxca,NA))
            szpft = NULL
            patch = with( data = datum$patch
                        , expr = mapply( FUN      = function(h,qw) ifelse(qw%!=%0,h/qw,NA)
                                       , h        = hflxca
                                       , qw       = qwflxca
                                       , SIMPLIFY = FALSE
                                       )#end mapply
                        )#end with
         }else if (this.vnam %in% "tratio"){
            emean = with(datum$emean,ifelse(wflxca %!=% 0,transp/wflxca,NA))
            szpft = NULL
            patch = with( data = datum$patch
                        , expr = mapply( FUN      = function(tp,w) ifelse(w%!=%0,tp/w,NA)
                                       , tp       = transp
                                       , w        = wflxca
                                       , SIMPLIFY = FALSE
                                       )#end mapply
                        )#end with
         }else{
            emean = datum$emean [[this.vnam]]
            szpft = datum$szpft [[this.vnam]]
            patch = datum$patch [[this.vnam]]
         }#end if (this.vnam %in% "tot.soil.c")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Scaling variables.                                                          #
         #---------------------------------------------------------------------------------#
         areaco  = datum$cohort$area
         areapa  = datum$patch$area
         agepa   = datum$patch$age
         loreypa = datum$patch$can.depth
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find out how to aggregate the running monthly mean.                        #
         #---------------------------------------------------------------------------------#
         if ( (! is.null(emean)) && (is.null(dim(emean)))){
            emean = matrix(emean,ncol=1)
         }#end if (! is.null(datum$emean))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Make the dimensions of emean and szpft match tomonth.                        #
         #---------------------------------------------------------------------------------#
         if (! is.null(emean)){
            xx.emean     = emean
            emean        = matrix( data     = NA
                                 , nrow     = nemean
                                 , ncol     = ncol(xx.emean)
                                 , dimnames = list( as.character(model$tomonth)
                                                  , dimnames(xx.emean)[[2]]
                                                  )#end list
                                 )#end matrix
            use          = datum$when %wr% model$tomonth
            eidx         = match(datum$when[use],model$tomonth)
            emean[eidx,] = xx.emean[use,]
         }#end if (! is.null(emean))
         if (! is.null(szpft)){
            xx.szpft      = szpft
            szpft         = array( data     = NA
                                 , dim      = c(nemean,dim(xx.szpft)[c(2,3)])
                                 , dimnames = list( as.character(model$tomonth)
                                                  , c(dbhkeys,"all")
                                                  , pft$key
                                                  )#end list
                                  )#end matrix
            use           = datum$when %wr% model$tomonth
            eidx          = match(datum$when[use],model$tomonth)
            szpft[eidx,,] = xx.szpft[use,,]
         }#end if (! is.null(szpft))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the annual mean, and crop the monthly mean to the period of interest.  #
         #---------------------------------------------------------------------------------#
         if (! is.null(emean)){
            ymean       = qapply( X     = emean
                                , INDEX = numyears(model$tomonth)
                                , DIM   = 1
                                , FUN   = mean
                                , na.rm = FALSE
                                )#end qapply
            mmean       = matrix(data=NA,nrow=12,ncol=ncol(emean))
            first       = qapply( X     = emean
                                , INDEX = nummonths(model$tomonth)
                                , DIM   = 1
                                , FUN   = mean
                                , na.rm = TRUE
                                )#end qapply
            idx         = as.numeric(dimnames(first)[[1]])
            mmean[idx,] = first
         }else{
            ymean = NULL
            mmean = NULL
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Crop the monthly means by PFT and DBH, keeping only the period of interest. #
         #---------------------------------------------------------------------------------#
         if (! is.null(szpft)){
            em.szpft        = szpft
            ym.szpft        = qapply( X     = em.szpft
                                    , INDEX = numyears(model$tomonth)
                                    , DIM   = 1
                                    , FUN   = mean
                                    , na.rm = FALSE
                                    )#end qapply
            mm.szpft        = array(data=NA,dim=c(12,dim(em.szpft)[c(2,3)]))
            mm.first        = qapply( X     = em.szpft
                                    , INDEX = nummonths(model$tomonth)
                                    , DIM   = 1
                                    , FUN   = mean
                                    , na.rm = TRUE
                                    )#end qapply
            idx             = as.numeric(dimnames(mm.first)[[1]])
            mm.szpft[idx,,] = mm.first

         }else{
            em.szpft = NULL
            ym.szpft = NULL
            mm.szpft = NULL
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find mortality rates, "interest" style.                                     #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% c("mort","dimort","ncbmort","firemort")){
            if (! is.null(emean   )) emean    = 100. * ( 1.0 - exp( - emean    ) )
            if (! is.null(ymean   )) ymean    = 100. * ( 1.0 - exp( - ymean    ) )
            if (! is.null(mmean   )) mmean    = 100. * ( 1.0 - exp( - mmean    ) )
            if (! is.null(szpft   )) szpft    = 100. * ( 1.0 - exp( - szpft    ) )
            if (! is.null(em.szpft)) em.szpft = 100. * ( 1.0 - exp( - em.szpft ) )
            if (! is.null(ym.szpft)) ym.szpft = 100. * ( 1.0 - exp( - ym.szpft ) )
            if (! is.null(mm.szpft)) mm.szpft = 100. * ( 1.0 - exp( - mm.szpft ) )
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Aggregate the patch-level data for the last cycle ("equilibrium").          #
         # em.age   -- Equivalent to emean                                                 #
         # mm.age   -- Equivalent to mmean.                                                #
         # pdf.area -- PDF weighted by area and number of occurrences during the period.   #
         #---------------------------------------------------------------------------------#
         if (! is.null(patch)){
            em.age     = array(data=NA,dim=c(nemean,n.age.at  ))
            em.lorey   = array(data=NA,dim=c(nemean,n.lorey.at))
            pdf.val    = mapply(FUN=numeric,length=rep(0,times=12),SIMPLIFY=FALSE)
            pdf.wgt    = mapply(FUN=numeric,length=rep(0,times=12),SIMPLIFY=FALSE)
            for (e in emean.loop){
               now       = model$tomonth[e]
               mm        = nummonths(now)
               yyyy      = numyears (now)
               stamp     = paste0("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm))
               vnow      = patch [[stamp]]
               agenow    = pmax(1/6,agepa[[stamp]])
               lnagenow  = log(agenow)
               loreynow  = loreypa[[stamp]]
               areanow   = areapa[[stamp]]
               oa        = order(lnagenow)
               ol        = order(loreynow)

               #----- Interpolate patch properties to fixed age classes. ------------------#
               if (length(vnow) > 1){
                  lorey.fun    = approxfun(x=loreynow[ol],y=vnow[ol])
                  lnage.fun    = approxfun(x=lnagenow[oa],y=vnow[oa])
                  lorey.int    = lorey.fun(v=lorey.at)
                  lnage.int    = lnage.fun(v=lnage.at)
                  em.lorey[e,] = ifelse(lorey.at %wr% range(loreynow),lorey.int,NA)
                  em.age  [e,] = ifelse(lnage.at %wr% range(lnagenow),lnage.int,NA)
               }else{
                  em.lorey[e,] = weighted.mean(x=vnow,w=areanow)
                  em.age  [e,] = weighted.mean(x=vnow,w=areanow)
               }#end if
               #---------------------------------------------------------------------------#


               #----- Append values and weights (area divided by number of times). --------#
               pdf.val[[mm]] = c(pdf.val[[mm]],vnow          )
               pdf.wgt[[mm]] = c(pdf.wgt[[mm]],areanow/nemean)
               #---------------------------------------------------------------------------#
            }#end for (e in emean.loop)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Find the monthly means.                                                  #
            #------------------------------------------------------------------------------#
            mm.age   = qapply( X     = em.age
                             , INDEX = nummonths(model$tomonth)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end qapply
            mm.age   = ifelse(is.finite(mm.age),mm.age,NA)
            mm.lorey = qapply( X     = em.lorey
                             , INDEX = nummonths(model$tomonth)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end qapply
            mm.lorey = ifelse(is.finite(mm.lorey),mm.lorey,NA)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Find the annual means.                                                   #
            #------------------------------------------------------------------------------#
            ym.age   = qapply( X     = em.age
                             , INDEX = numyears(model$tomonth)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = FALSE
                             )#end qapply
            ym.age   = ifelse(is.finite(ym.age),ym.age,NA)
            ym.lorey = qapply( X     = em.lorey
                             , INDEX = numyears(model$tomonth)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = FALSE
                             )#end qapply
            ym.lorey = ifelse(is.finite(ym.lorey),ym.lorey,NA)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Find the seasonal means.                                                 #
            #------------------------------------------------------------------------------#
            zm.age   = qapply( X     = em.age
                             , INDEX = season(model$tomonth,add.year=FALSE)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end qapply
            zm.age   = ifelse(is.finite(zm.age),zm.age,NA)
            zm.lorey = qapply( X     = em.lorey
                             , INDEX = season(model$tomonth,add.year=FALSE)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end qapply
            zm.lorey = ifelse(is.finite(zm.lorey),zm.lorey,NA)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Find the density function for each month.                                #
            #------------------------------------------------------------------------------#
            mm.pdf = t(mapply( FUN      = function(x,weights=NULL,dens.min=0.0001,...){
                                             out = density.safe(x,weights,...)$y
                                             return(out)
                                          }#end function
                             , x        = pdf.val
                             , weights  = pdf.wgt
                             , MoreArgs = list( n    = n.dens
                                              , from = min(var.range[1,v,,],na.rm=TRUE)
                                              , to   = max(var.range[2,v,,],na.rm=TRUE)
                                              )#end list
                             , SIMPLIFY = TRUE
                             )#end mapply
                      )#end t
            #------------------------------------------------------------------------------#
         }else{
            em.age   = NULL
            mm.age   = NULL
            ym.age   = NULL
            zm.age   = NULL
            em.lorey = NULL
            mm.lorey = NULL
            ym.lorey = NULL
            zm.lorey = NULL
            mm.pdf   = NULL
         }#end if (! is.null(patch))
         #---------------------------------------------------------------------------------#


         #----- Save variables to a list. -------------------------------------------------#
         model[[this.vnam]] = list( emean    = emean
                                  , ymean    = ymean
                                  , mmean    = mmean
                                  , em.szpft = em.szpft
                                  , ym.szpft = ym.szpft
                                  , mm.szpft = mm.szpft
                                  , em.age   = em.age
                                  , mm.age   = mm.age
                                  , ym.age   = ym.age
                                  , zm.age   = zm.age
                                  , em.lorey = em.lorey
                                  , mm.lorey = mm.lorey
                                  , ym.lorey = ym.lorey
                                  , zm.lorey = zm.lorey
                                  , mm.pdf   = mm.pdf
                                  )#end list
         rm(list=c("emean","ymean","mmean","szpft","em.szpft","ym.szpft","mm.szpft"
                  ,"em.age","mm.age","zm.age","ym.age"
                  ,"em.lorey","mm.lorey","zm.lorey","ym.lorey","mm.pdf")
           )#end rm
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ncompvar))
      #------------------------------------------------------------------------------------#



      #----- Save the data and free some memory. ------------------------------------------#
      this.sim[[simul$name[s]]] = model
      rm(list=c("model"))
      #------------------------------------------------------------------------------------#
   }#end for (s in sequence(nsimul))
   #---------------------------------------------------------------------------------------#





   #----- Copy the data to the results. ---------------------------------------------------#
   res.iata    = paste("res",iata,sep=".")
   assign(res.iata,this.sim)
   res[[iata]] = this.sim
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      Save processed data to RData.                                                    #
   #---------------------------------------------------------------------------------------#
   rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
   cat0(" + Saving processed data to ",basename(rdata.iata),"...")
   dummy = save(list=c(res.iata), file=rdata.iata)
   rm(res.iata,this.sim)
   #---------------------------------------------------------------------------------------#

}#end for (p in sequence(loop.nsites))
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
#     Plot the monthly mean time series for each site.                                     #
#------------------------------------------------------------------------------------------#
if (plot.ts.emean){
   cat0(" + Plot time series of monthly means...")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      #------------------------------------------------------------------------------------#



      cat0("   - ",this.desc,"...")
      #------------------------------------------------------------------------------------#
      #     Loop over all sites and simulations, and get the range.                        #
      #------------------------------------------------------------------------------------#
      em.xrange = array(NA,dim=c(2,nsites,nsimul))
      em.yrange = array(NA,dim=c(2,nsites,nsimul))
      for (p in sequence(nsites)){
         iata = sites$iata[p]
         for (s in sequence(nsimul)){
            #------ Get the data. ---------------------------------------------------------#
            sname    = simul$name[s]
            model    = res[[iata]][[sname]]
            tomonth  = model$tomonth
            toyear   = model$toyear
            emean    = model[[this.vnam]]$emean
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Update range.                                                            #
            #------------------------------------------------------------------------------#
            em.xrange[,p,s] = range(c(em.xrange[,p,s],tomonth),finite=TRUE)
            em.yrange[,p,s] = range(c(em.yrange[,p,s],emean  ),finite=TRUE)
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Make one panel for each simulation, and one plot per site.                     #
      #------------------------------------------------------------------------------------#
      cat0("     * Plot by sites...")
      for (p in sequence(nsites)){
         #----- Get the basic information. ------------------------------------------------#
         iata            = sites$iata[p]
         this.longname   = sites$desc[p]
         em.xlimit       = pretty.xylim(u=em.xrange[,p,],fracexp=0.0,is.log=FALSE)
         em.ylimit       = pretty.xylim(u=em.yrange[,p,],fracexp=0.0,is.log=FALSE)
         em.pretty       = pretty.time(chron(em.xlimit),n=5)
         cat0("       > ",this.longname,"...")
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = this.longname
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         lex     = "Time"
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ts.emean
            fichier = file.path( out.now
                               , paste0("ts_emean-",this.vnam,"-",iata,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=ssize$width,height=ssize$height
                         ,pointsize=ptsz,paper=ssize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                  ,pointsize=ptsz,paper=ssize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            layout(mat= rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = simleg.key
                  , fill    = simcol.key
                  , border  = simcol.key
                  , xpd     = TRUE
                  , bty     = "n"
                  , cex     = 0.7
                  )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot monthly means.                                                     #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,3.1,1.6))
            plot.new()
            plot.window(xlim=em.xlimit,ylim=em.ylimit)
            axis(side=1,las=1,at=em.pretty$levels,labels=em.pretty$labels)
            axis(side=2,las=2)
            title(main=letitre,xlab=lex,ylab=ley,cex.main=1.0)
            for (s in sequence(nsimul)){
               #----- Load variables. -----------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               tomonth  = model$tomonth
               emean    = model[[this.vnam]]$emean
               lines(x=tomonth,y=emean,col=simcol.key[s],lwd=1.5,type="l")
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
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.emean)
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
#     Plot the time series of annual means.                                                #
#------------------------------------------------------------------------------------------#
if (plot.ts.ymean){
   cat0(" + Plot time series of annual means...")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      #------------------------------------------------------------------------------------#



      cat0("   - ",this.desc,"...")
      #------------------------------------------------------------------------------------#
      #     Loop over all sites and simulations, and get the range.                        #
      #------------------------------------------------------------------------------------#
      ym.xrange = array(NA,dim=c(2,nsites,nsimul))
      ym.yrange = array(NA,dim=c(2,nsites,nsimul))
      for (p in sequence(nsites)){
         iata = sites$iata[p]
         for (s in sequence(nsimul)){
            #------ Get the data. ---------------------------------------------------------#
            sname    = simul$name[s]
            model    = res[[iata]][[sname]]
            toyear   = model$toyear
            ymean    = model[[this.vnam]]$ymean
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Update range.                                                            #
            #------------------------------------------------------------------------------#
            ym.xrange[,p,s] = range(c(ym.xrange[,p,s],toyear ),finite=TRUE)
            ym.yrange[,p,s] = range(c(ym.yrange[,p,s],ymean  ),finite=TRUE)
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Make one panel for each simulation, and one plot per site.                     #
      #------------------------------------------------------------------------------------#
      cat0("     * Plot by sites...")
      for (p in sequence(nsites)){
         #----- Get the basic information. ------------------------------------------------#
         iata            = sites$iata[p]
         this.longname   = sites$desc[p]
         ym.xlimit       = pretty.xylim(u=ym.xrange[,p,],fracexp=0.0,is.log=FALSE)
         ym.ylimit       = pretty.xylim(u=ym.yrange[,p,],fracexp=0.0,is.log=FALSE)
         cat0("       > ",this.longname,"...")
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = this.longname
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         lex     = "Time"
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ts.ymean
            fichier = file.path( out.now
                               , paste0("ts_ymean-",this.vnam,"-",iata,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=ssize$width,height=ssize$height
                         ,pointsize=ptsz,paper=ssize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                  ,pointsize=ptsz,paper=ssize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            layout(mat= rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = simleg.key
                  , fill    = simcol.key
                  , border  = simcol.key
                  , xpd     = TRUE
                  , bty     = "n"
                  , cex     = 0.7
                  )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot annual means.                                                      #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,3.1,1.6))
            plot.new()
            plot.window(xlim=ym.xlimit,ylim=ym.ylimit)
            axis(side=1,las=1)
            axis(side=2,las=2)
            title(main="Annual means",xlab=lex,ylab=ley,cex.main=1.0)
            for (s in sequence(nsimul)){
               #----- Load variables. -----------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               toyear   = model$toyear
               ymean    = model[[this.vnam]]$ymean
               lines(x=toyear,y=ymean,col=simcol.key[s],lwd=1.5,type="l")
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
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.ymean)
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
#     Plot the time series by PFT for each site.                                           #
#------------------------------------------------------------------------------------------#
if (plot.ts.pft){
   cat0(" + Plot long-term time series by PFT...")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.szpft        = this.compvar$szpftvar
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.szpft){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         xrange = array(NA,dim=c(2,nsites,nsimul,npft+1))
         yrange = array(NA,dim=c(2,nsites,nsimul,npft+1))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               toyear   = model$toyear
               ym.szpft = model[[this.vnam]]$ym.szpft
               nsize    = dim(ym.szpft)[2]
               
               for (f in pft.use){
                  #------------------------------------------------------------------------#
                  #     Update range.                                                      #
                  #------------------------------------------------------------------------#
                  ym.pft         = ym.szpft[,nsize,f]
                  xrange[,p,s,f] = range(c(xrange[,p,s,f],toyear),finite=TRUE)
                  yrange[,p,s,f] = range(c(yrange[,p,s,f],ym.pft),finite=TRUE)
                  #------------------------------------------------------------------------#
               }#end for (f in pft.use)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#




            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"Annual means",sep=" - ")
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            lex     = "Year"
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ts.pft$variables
               fichier = file.path( out.now
                                  , paste0("ts_pft-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=fsize$width,height=fsize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=fsize$width,height=fsize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=fsize$width,height=fsize$height
                            ,pointsize=ptsz,paper=fsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=fsize$width,height=fsize$height
                     ,pointsize=ptsz,paper=fsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout( mat     = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                     , heights = c(rep((1.-f.leg)/lo.pft$nrow,lo.pft$nrow),f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#


               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = simleg.key
                     , fill    = simcol.key
                     , border  = simcol.key
                     , ncol    = min(3,pretty.box(n=n.sim)$ncol)
                     , title   = expression(bold("Simulation"))
                     , xpd     = TRUE
                     , bty     = "o"
                     , cex     = 0.8
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over PFTs then by simulations.                                  #
               #---------------------------------------------------------------------------#
               for (f in pft.use){
                  xlimit = pretty.xylim(u=xrange[,p,,f],fracexp=0.0,is.log=FALSE)
                  ylimit = pretty.xylim(u=yrange[,p,,f],fracexp=0.0,is.log=FALSE)

                  #----- Open window and plot all time series by simulation. --------------#
                  par(mar=c(4.1,4.6,2.1,0.6))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  for (s in sequence(nsimul)){
                     #----- Load variables. -----------------------------------------------#
                     sname    = simul$name[s]
                     model    = res[[iata]][[sname]]
                     toyear   = model$toyear
                     ym.szpft = model[[this.vnam]]$ym.szpft
                     nsize    = dim(ym.szpft)[2]
                     #---------------------------------------------------------------------#

                     lines( x    = toyear
                          , y    = ym.szpft[,nsize,f]
                          , col  = simcol.key[s]
                          , lwd  = 2.0
                          , type = "l"
                          )#end lines
                  }#end for (s in sequence(nsimul))
                  axis(side=1,las=1)
                  axis(side=2,las=1)
                  title(main=pft$name[f],line=0.5)
                  box()
                  #------------------------------------------------------------------------#
               }#end for (f in pft.use)
               #---------------------------------------------------------------------------#



               #----- Plot the global title. ----------------------------------------------#
               gtitle( main      = letitre
                     , xlab      = lex
                     , ylab      = ley
                     , off.xlab  = f.leg / (1 + f.leg)
                     , cex.main  = 1.0
                     , line.main = 2.5
                     )#end gtitle
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot default simulation.                                                    #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plotting default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname           = simul$name[s]
         sdesc           = simul$desc[s]
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = paste(sdesc,"Annual means",sep=" - ")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         lex     = "Year"
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ts.pft$default
            fichier = file.path( out.now
                               , paste0("ts_pft-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=fsize$width,height=fsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=fsize$width,height=fsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=fsize$width,height=fsize$height
                         ,pointsize=ptsz,paper=fsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=fsize$width,height=fsize$height
                  ,pointsize=ptsz,paper=fsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0.25,2.5,0))
            layout( mat     = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                  , heights = c(rep((1.-f.leg)/lo.pft$nrow,lo.pft$nrow),f.leg)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = sites$iata
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(n=nsites)$ncol)
                  , title   = expression(bold("Sites"))
                  , xpd     = TRUE
                  , bty     = "o"
                  , cex     = 0.8
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over PFTs.                                                         #
            #------------------------------------------------------------------------------#
            for (f in pft.use){
               xlimit = pretty.xylim(u=xrange[,,s,f],fracexp=0.0,is.log=FALSE)
               ylimit = pretty.xylim(u=yrange[,,s,f],fracexp=0.0,is.log=FALSE)


               #----- Open window and plot all time series for this PFT. ------------------#
               par(mar=c(4.1,4.6,2.1,0.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     Plot all sites.                                                       #
               #---------------------------------------------------------------------------#
               for (p in sequence(nsites)){
                  iata            = sites$iata[p]
                  this.longname   = sites$desc[p]
                  model           = res[[iata]][[sname]]
                  toyear          = model$toyear
                  ym.szpft        = model[[this.vnam]]$ym.szpft
                  nsize           = dim(ym.szpft)[2]
                  lines( x    = toyear
                       , y    = ym.szpft[,nsize,f]
                       , col  = sites$col[p]
                       , lwd  = 2.0
                       , type = "l"
                       )#end lines
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot axis-related stuff. --------------------------------------------#
               axis(side=1,las=1)
               axis(side=2,las=1)
               title(main=pft$name[f],line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (f in pft.use)
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , xlab      = lex
                  , ylab      = ley
                  , off.xlab  = f.leg / (1. + f.leg)
                  , line.main = 2.5
                  , line.ylab = 3.0
                  , cex.axis  = 1.0
                  )#end gtitle
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
      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.pft)
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
#     Plot the time series by DBH for each site.                                           #
#------------------------------------------------------------------------------------------#
if (plot.ts.dbh){
   cat0(" + Plot long-term time series by DBH...")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.szpft        = this.compvar$szpftvar
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.szpft){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         xrange = array(NA,dim=c(2,nsites,nsimul,ndbh))
         yrange = array(NA,dim=c(2,nsites,nsimul,ndbh))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               toyear   = model$toyear
               ym.szpft = model[[this.vnam]]$ym.szpft
               nsize    = dim(ym.szpft)[3]

               for (f in dbh.use){
                  #------------------------------------------------------------------------#
                  #     Update range.                                                      #
                  #------------------------------------------------------------------------#
                  ym.dbh         = ym.szpft[,f,nsize]
                  xrange[,p,s,f] = range(c(xrange[,p,s,f],toyear),finite=TRUE)
                  yrange[,p,s,f] = range(c(yrange[,p,s,f],ym.dbh),finite=TRUE)
                  #------------------------------------------------------------------------#
               }#end for (f in dbh.use)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#




            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"Annual means",sep=" - ")
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            lex     = "Year"
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ts.dbh$variables
               fichier = file.path( out.now
                                  , paste0("ts_dbh-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=dsize$width,height=dsize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=dsize$width,height=dsize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=dsize$width,height=dsize$height
                            ,pointsize=ptsz,paper=dsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=dsize$width,height=dsize$height
                     ,pointsize=ptsz,paper=dsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout( mat     = rbind(lo.dbh$mat.off,rep(1,times=lo.dbh$ncol))
                     , heights = c(rep((1.-f.leg)/lo.dbh$nrow,lo.dbh$nrow),f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#


               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = simleg.key
                     , fill    = simcol.key
                     , border  = simcol.key
                     , ncol    = min(3,pretty.box(n=n.sim)$ncol)
                     , title   = expression(bold("Simulation"))
                     , xpd     = TRUE
                     , bty     = "o"
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over DBHs then by simulations.                                  #
               #---------------------------------------------------------------------------#
               for (f in dbh.use){
                  xlimit = pretty.xylim(u=xrange[,p,,f],fracexp=0.0,is.log=FALSE)
                  ylimit = pretty.xylim(u=yrange[,p,,f],fracexp=0.0,is.log=FALSE)

                  #----- Open window and plot all time series by simulation. --------------#
                  par(mar=c(4.1,4.6,2.1,0.6))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  for (s in sequence(nsimul)){
                     #----- Load variables. -----------------------------------------------#
                     sname    = simul$name[s]
                     model    = res[[iata]][[sname]]
                     toyear   = model$toyear
                     ym.szpft = model[[this.vnam]]$ym.szpft
                     nsize    = dim(ym.szpft)[3]
                     #---------------------------------------------------------------------#

                     lines( x    = toyear
                          , y    = ym.szpft[,f,nsize]
                          , col  = simcol.key[s]
                          , lwd  = 2.0
                          , type = "l"
                          )#end lines
                  }#end for (s in sequence(nsimul))
                  axis(side=1,las=1)
                  axis(side=2,las=1)
                  title(main=dbhnames[f],line=0.5)
                  box()
                  #------------------------------------------------------------------------#
               }#end for (f in pft.use)
               #---------------------------------------------------------------------------#



               #----- Plot the global title. ----------------------------------------------#
               gtitle( main      = letitre
                     , xlab      = lex
                     , ylab      = ley
                     , off.xlab  = f.leg / (1 + f.leg)
                     , cex.main  = 1.0
                     , line.main = 2.5
                     )#end gtitle
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot default simulation.                                                    #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plotting default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname           = simul$name[s]
         sdesc           = simul$desc[s]
         xlimit          = pretty.xylim(u=xrange[,,s,f],fracexp=0.0,is.log=FALSE)
         ylimit          = pretty.xylim(u=yrange[,,s,f],fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = paste(sdesc,"Annual means",sep=" - ")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         lex     = "Year"
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ts.dbh$default
            fichier = file.path( out.now
                               , paste0("ts_dbh-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=dsize$width,height=dsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=dsize$width,height=dsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=dsize$width,height=dsize$height
                         ,pointsize=ptsz,paper=dsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=dsize$width,height=dsize$height
                  ,pointsize=ptsz,paper=dsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0.25,2.5,0))
            layout( mat     = rbind(lo.dbh$mat.off,rep(1,times=lo.dbh$ncol))
                  , heights = c(rep((1.-f.leg)/lo.dbh$nrow,lo.dbh$nrow),f.leg)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = sites$iata
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(n=nsites)$ncol)
                  , title   = expression(bold("Sites"))
                  , xpd     = TRUE
                  , bty     = "o"
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over DBH classes.                                                  #
            #------------------------------------------------------------------------------#
            for (f in dbh.use){
               xlimit = pretty.xylim(u=xrange[,,s,f],fracexp=0.0,is.log=FALSE)
               ylimit = pretty.xylim(u=yrange[,,s,f],fracexp=0.0,is.log=FALSE)


               #----- Open window and plot all time series for this PFT. ------------------#
               par(mar=c(4.1,4.6,2.1,0.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot all sites.                                                       #
               #---------------------------------------------------------------------------#
               for (p in sequence(nsites)){
                  iata            = sites$iata[p]
                  this.longname   = sites$desc[p]
                  model           = res[[iata]][[sname]]
                  toyear          = model$toyear
                  ym.szpft        = model[[this.vnam]]$ym.szpft
                  nsize           = dim(ym.szpft)[3]
                  lines( x    = toyear
                       , y    = ym.szpft[,f,nsize]
                       , col  = sites$col[p]
                       , lwd  = 2.0
                       , type = "l"
                       )#end lines
               }#end for (s in sequence(nsimul))
               axis(side=1,las=1)
               axis(side=2,las=1)
               title(main=dbhnames[f],line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (f in pft.use)
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , xlab      = lex
                  , ylab      = ley
                  , off.xlab  = f.leg / (1. + f.leg)
                  , line.main = 2.5
                  , line.ylab = 3.0
                  , cex.axis  = 1.0
                  )#end gtitle
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



      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.pft)
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
#     Plot the monthly means by PFT for each site.                                         #
#------------------------------------------------------------------------------------------#
if (plot.mm.pft){
   cat0(" + Plot monthly means by PFT...")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.szpft        = this.compvar$szpftvar
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.szpft){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         yrange = array(NA,dim=c(2,nsites,nsimul,npft+1))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               mm.szpft = model[[this.vnam]]$mm.szpft
               nsize    = dim(mm.szpft)[2]
               for (f in pft.use){
                  #------------------------------------------------------------------------#
                  #     Update range.                                                      #
                  #------------------------------------------------------------------------#
                  mm.pft         = mm.szpft[,nsize,f]
                  yrange[,p,s,f] = range(c(yrange[,p,s,f],mm.pft),finite=TRUE)
                  #------------------------------------------------------------------------#
               }#end for (f in pft.use)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Month limits.                                                               #
         #---------------------------------------------------------------------------------#
         xlimit          = c(0.5,12.5)
         xat             = seq_along(month.abb)
         xlabels         = substring(month.abb,1,1)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#




            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"Monthly means",sep=" - ")
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$mm.pft$variables
               fichier = file.path( out.now
                                  , paste0("mm_pft-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=fsize$width,height=fsize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=fsize$width,height=fsize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=fsize$width,height=fsize$height
                            ,pointsize=ptsz,paper=fsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=fsize$width,height=fsize$height
                     ,pointsize=ptsz,paper=fsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout( mat     = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                     , heights = c(rep((1.-f.leg)/lo.pft$nrow,lo.pft$nrow),f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#


               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = simleg.key
                     , fill    = simcol.key
                     , border  = simcol.key
                     , ncol    = min(3,pretty.box(n=n.sim)$ncol)
                     , title   = expression(bold("Simulation"))
                     , xpd     = TRUE
                     , bty     = "o"
                     , cex     = 0.8
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over PFTs then by simulations.                                  #
               #---------------------------------------------------------------------------#
               for (f in pft.use){
                  ylimit = pretty.xylim(u=yrange[,p,,f],fracexp=0.0,is.log=FALSE)

                  #----- Open window and plot all time series by simulation. --------------#
                  par(mar=c(4.1,4.6,2.1,0.6))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  for (s in sequence(nsimul)){
                     #----- Load variables. -----------------------------------------------#
                     sname    = simul$name[s]
                     model    = res[[iata]][[sname]]
                     toyear   = model$toyear
                     mm.szpft = model[[this.vnam]]$mm.szpft
                     nsize    = dim(mm.szpft)[2]
                     #---------------------------------------------------------------------#

                     lines( x    = xat
                          , y    = mm.szpft[,nsize,f]
                          , col  = simcol.key[s]
                          , lwd  = 2.0
                          , type = "l"
                          )#end lines
                  }#end for (s in sequence(nsimul))
                  axis(side=1,las=1,at=xat,labels=xlabels)
                  axis(side=2,las=1)
                  title(main=pft$name[f],line=0.5)
                  box()
                  #------------------------------------------------------------------------#
               }#end for (f in pft.use)
               #---------------------------------------------------------------------------#



               #----- Plot the global title. ----------------------------------------------#
               gtitle( main      = letitre
                     , ylab      = ley
                     , cex.main  = 1.0
                     , line.main = 2.5
                     )#end gtitle
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot default simulation.                                                    #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plotting default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname           = simul$name[s]
         sdesc           = simul$desc[s]
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = paste(sdesc,"Monthly means",sep=" - ")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$mm.pft$default
            fichier = file.path( out.now
                               , paste0("mm_pft-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=fsize$width,height=fsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=fsize$width,height=fsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=fsize$width*depth,height=fsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=fsize$width,height=fsize$height
                         ,pointsize=ptsz,paper=fsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=fsize$width,height=fsize$height
                  ,pointsize=ptsz,paper=fsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0.25,2.5,0))
            layout( mat     = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                  , heights = c(rep((1.-f.leg)/lo.pft$nrow,lo.pft$nrow),f.leg)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = sites$iata
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(n=nsites)$ncol)
                  , title   = expression(bold("Sites"))
                  , xpd     = TRUE
                  , bty     = "o"
                  , cex     = 0.8
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over PFTs.                                                         #
            #------------------------------------------------------------------------------#
            for (f in pft.use){
               ylimit = pretty.xylim(u=yrange[,,s,f],fracexp=0.0,is.log=FALSE)


               #----- Open window and plot all time series for this PFT. ------------------#
               par(mar=c(4.1,4.6,2.1,0.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot all sites.                                                       #
               #---------------------------------------------------------------------------#
               for (p in sequence(nsites)){
                  iata            = sites$iata[p]
                  this.longname   = sites$desc[p]
                  model           = res[[iata]][[sname]]
                  mm.szpft        = model[[this.vnam]]$mm.szpft
                  nsize           = dim(mm.szpft)[2]
                  lines( x    = xat
                       , y    = mm.szpft[,nsize,f]
                       , col  = sites$col[p]
                       , lwd  = 2.0
                       , type = "l"
                       )#end lines
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot axis-related stuff. --------------------------------------------#
               axis(side=1,las=1,at=xat,labels=xlabels)
               axis(side=2,las=1)
               title(main=pft$name[f],line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (f in pft.use)
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , ylab      = ley
                  , line.main = 2.5
                  , line.ylab = 3.0
                  , cex.axis  = 1.0
                  )#end gtitle
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
      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.pft)
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
#     Plot monthly means by DBH for each site.                                             #
#------------------------------------------------------------------------------------------#
if (plot.mm.dbh){
   cat0(" + Plot monthly means by DBH...")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.szpft        = this.compvar$szpftvar
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.szpft){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         yrange = array(NA,dim=c(2,nsites,nsimul,ndbh))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               toyear   = model$toyear
               mm.szpft = model[[this.vnam]]$mm.szpft
               nsize    = dim(mm.szpft)[3]
               #---------------------------------------------------------------------------#
               for (f in dbh.use){
                  #------------------------------------------------------------------------#
                  #     Update range.                                                      #
                  #------------------------------------------------------------------------#
                  mm.dbh         = mm.szpft[,f,nsize]
                  yrange[,p,s,f] = range(c(yrange[,p,s,f],mm.dbh),finite=TRUE)
                  #------------------------------------------------------------------------#
               }#end for (f in dbh.use)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Month limits.                                                               #
         #---------------------------------------------------------------------------------#
         xlimit          = c(0.5,12.5)
         xat             = seq_along(month.abb)
         xlabels         = substring(month.abb,1,1)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#




            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"Monthly means",sep=" - ")
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$mm.dbh$variables
               fichier = file.path( out.now
                                  , paste0("mm_dbh-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=dsize$width,height=dsize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=dsize$width,height=dsize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=dsize$width,height=dsize$height
                            ,pointsize=ptsz,paper=dsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=dsize$width,height=dsize$height
                     ,pointsize=ptsz,paper=dsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout( mat     = rbind(lo.dbh$mat.off,rep(1,times=lo.dbh$ncol))
                     , heights = c(rep((1.-f.leg)/lo.dbh$nrow,lo.dbh$nrow),f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#


               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = simleg.key
                     , fill    = simcol.key
                     , border  = simcol.key
                     , ncol    = min(3,pretty.box(n=n.sim)$ncol)
                     , title   = expression(bold("Simulation"))
                     , xpd     = TRUE
                     , bty     = "o"
                     , cex     = 0.8
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over DBHs then by simulations.                                  #
               #---------------------------------------------------------------------------#
               for (f in dbh.use){
                  ylimit = pretty.xylim(u=yrange[,p,,f],fracexp=0.0,is.log=FALSE)

                  #----- Open window and plot all time series by simulation. --------------#
                  par(mar=c(4.1,4.6,2.1,0.6))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  for (s in sequence(nsimul)){
                     #----- Load variables. -----------------------------------------------#
                     sname    = simul$name[s]
                     model    = res[[iata]][[sname]]
                     mm.szpft = model[[this.vnam]]$mm.szpft
                     nsize    = dim(mm.szpft)[3]
                     #---------------------------------------------------------------------#

                     lines( x    = xat
                          , y    = mm.szpft[,f,nsize]
                          , col  = simcol.key[s]
                          , lwd  = 2.0
                          , type = "l"
                          )#end lines
                  }#end for (s in sequence(nsimul))
                  axis(side=1,las=1,at=xat,labels=xlabels)
                  axis(side=2,las=1)
                  title(main=dbhnames[f],line=0.5)
                  box()
                  #------------------------------------------------------------------------#
               }#end for (f in pft.use)
               #---------------------------------------------------------------------------#



               #----- Plot the global title. ----------------------------------------------#
               gtitle( main      = letitre
                     , ylab      = ley
                     , cex.main  = 1.0
                     , line.main = 2.5
                     )#end gtitle
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot default simulation.                                                    #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plotting default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname           = simul$name[s]
         sdesc           = simul$desc[s]
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = paste(sdesc,"Monthly means",sep=" - ")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$mm.dbh$default
            fichier = file.path( out.now
                               , paste0("mm_dbh-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=dsize$width,height=dsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=dsize$width,height=dsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=dsize$width*depth,height=dsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=dsize$width,height=dsize$height
                         ,pointsize=ptsz,paper=dsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=dsize$width,height=dsize$height
                  ,pointsize=ptsz,paper=dsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0.25,2.5,0))
            layout( mat     = rbind(lo.dbh$mat.off,rep(1,times=lo.dbh$ncol))
                  , heights = c(rep((1.-f.leg)/lo.dbh$nrow,lo.dbh$nrow),f.leg)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = sites$iata
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(n=length(nsites))$ncol)
                  , title   = expression(bold("Sites"))
                  , xpd     = TRUE
                  , bty     = "o"
                  , cex     = 0.8
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over DBH classes.                                                  #
            #------------------------------------------------------------------------------#
            for (f in dbh.use){
               ylimit = pretty.xylim(u=yrange[,,s,f],fracexp=0.0,is.log=FALSE)


               #----- Open window and plot all time series for this PFT. ------------------#
               par(mar=c(4.1,4.6,2.1,0.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot all sites.                                                       #
               #---------------------------------------------------------------------------#
               for (p in sequence(nsites)){
                  iata            = sites$iata[p]
                  this.longname   = sites$desc[p]
                  model           = res[[iata]][[sname]]
                  mm.szpft        = model[[this.vnam]]$mm.szpft
                  nsize           = dim(mm.szpft)[3]
                  lines( x    = xat
                       , y    = mm.szpft[,f,nsize]
                       , col  = sites$col[p]
                       , lwd  = 2.0
                       , type = "l"
                       )#end lines
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot axis-related stuff. --------------------------------------------#
               axis(side=1,las=1,at=xat,labels=xlabels)
               axis(side=2,las=1)
               title(main=dbhnames[f],line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (f in dbh.use)
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , ylab      = ley
                  , line.main = 2.5
                  , line.ylab = 3.0
                  , cex.axis  = 1.0
                  )#end gtitle
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



      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.mm.pft)
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
#     Plot the mean annual cycle by age structure.                                         #
#------------------------------------------------------------------------------------------#
if (plot.ym.patch){
   cat0(" + Plot average properties as a function of age...")


   if (patch.aggr %in% "lorey"){
      x.at     = pretty(x=lorey.at)
      x.labels = sprintf("%g",x.at)
      xlimit   = range(lorey.at)
   }else{
      x.at     = pretty.log(x=age.at)
      x.labels = sprintf("%g",x.at)
      xlimit   = range(age.at)
   }#end if

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.patch        = this.compvar$patchvar
      zlog            = this.compvar$zlog
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.patch){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         yrange = array(NA,dim=c(2,nsites,nsimul))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               if (patch.aggr %in% "lorey"){
                  ym.patch   = colMeans(model[[this.vnam]]$ym.lorey)
               }else{
                  ym.patch   = colMeans(model[[this.vnam]]$ym.age)
               }#end if
               ym.patch = ifelse(is.finite(ym.patch),ym.patch,NA)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               if (zlog){
                  ym.patch     = ifelse(ym.patch %>% 0.0,ym.patch,NA)
                  yrange[,p,s] = range(c(yrange[,p,s],ym.patch),finite=TRUE)
               }else if (this.vnam %in% "bowen"){
                  ym.patch     = pmax(bmn,pmin(bmx,ym.patch))
                  yrange[,p,s] = range(c(yrange[,p,s],ym.patch),finite=TRUE)
               }else if (this.vnam %in% "tratio"){
                  ym.patch     = pmax(0.0,pmin(1.0,ym.patch))
                  yrange[,p,s] = range(c(yrange[,p,s],ym.patch),finite=TRUE)
               }else{
                  yrange[,p,s] = range(c(yrange[,p,s],ym.patch),finite=TRUE)
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #------ Make the scale for budget variables symmetric. ---------------------------#
         if (this.vnam %in% c("nee","nep","cba")){
            yrange = apply( X      = yrange
                          , MARGIN = c(2,3)
                          , FUN    = function(x) c(-1,1)*max(abs(x),na.rm=TRUE)
                          )#end apply
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata           = sites$iata[p]
            this.longname  = sites$desc[p]
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Grab data and arrange them by list.                                     #
            #------------------------------------------------------------------------------#
            xdat = list()
            ydat = list()
            for (s in sequence(nsimul)){
               sname        = simul$name[s]
               model        = res[[iata]][[sname]]
               if (patch.aggr %in% "lorey"){
                  xdat[[s]] = lorey.at
                  ym.patch  = colMeans(model[[this.vnam]]$ym.lorey)
               }else{
                  xdat[[s]] = age.at
                  ym.patch  = colMeans(model[[this.vnam]]$ym.age)
               }#end if
               ym.patch = ifelse(is.finite(ym.patch),ym.patch,NA)
               #---------------------------------------------------------------------------#
               if (this.vnam %in% "bowen"){
                  ydat[[s]] = pmax(bmn,pmin(bmx,ym.patch))
               }else if (this.vnam %in% "tratio"){
                  ydat[[s]] = pmax(0.0,pmin(1.0,ym.patch))
               }else{
                  ydat[[s]] = ym.patch
               }#end if
            }#end for (s in sequence(nsimul))
            ylimit = pretty.xylim(u=yrange[,p,],fracexp=0.0,is.log=zlog)
            if (zlog){
               plog     = "y"
               y.at     = pretty.log(ylimit)
               y.labels = sprintf("%g",y.at)
            }else{
               plog     = ""
               y.at     = pretty(ylimit)
               y.labels = sprintf("%g",y.at)
            }#end if
            #------------------------------------------------------------------------------#


            #------ Set some common features. ---------------------------------------------#
            letitre = paste0(this.desc," - ",this.longname)
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            if (patch.aggr %in% "lorey"){
               lex  = desc.unit(desc="Lorey's height",unit=untab$m)
            }else{
               lex  = desc.unit(desc="Age",unit=untab$yr)
               plog = paste0("x",plog)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ym.patch$variables
               fichier = file.path( out.now
                                  , paste0("ym_patch-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
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



               #---------------------------------------------------------------------------#
               #      Split the device area into two.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot legend.                                                         #
               #---------------------------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = simul$desc
                     , fill    = simul$colour
                     , border  = simul$colour
                     , cex     = 0.75
                     , xpd     = TRUE
                     , bty     = "n"
                     )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Plot simulations.                                                    #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.6,3.1,1.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=plog)
               for (s in sequence(nsimul)){
                  lines(x=xdat[[s]],y=ydat[[s]],lwd=1.5,col=simul$colour[s])
               }#end for
               box()
               axis(side=1,at=x.at,labels=x.labels)
               axis(side=2,at=y.at,labels=y.labels,las=1)
               title(main=letitre,xlab=lex,ylab=ley)
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plot default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname     = simul$name[s]
         sdesc     = simul$desc[s]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Grab data and arrange them by list.                                        #
         #---------------------------------------------------------------------------------#
         xdat = list()
         ydat = list()
         for (p in sequence(nsites)){
            iata         = sites$iata[p]
            longname     = sites$desc[p]
            model        = res[[iata]][[sname]]
            if (patch.aggr %in% "lorey"){
               xdat[[p]] = lorey.at
               ym.patch  = colMeans(model[[this.vnam]]$ym.lorey)
            }else{
               xdat[[p]] = age.at
               ym.patch  = colMeans(model[[this.vnam]]$ym.age)
            }#end if (patch.aggr %in% "lorey")
            if (this.vnam %in% "bowen"){
               ydat[[p]] = pmax(bmn,pmin(bmx,ym.patch))
            }else if (this.vnam %in% "tratio"){
               ydat[[p]] = pmax(0.0,pmin(1.0,ym.patch))
            }else{
               ydat[[p]] = ym.patch
            }#end if
            ym.patch = ifelse(is.finite(ym.patch),ym.patch,NA)
         }#end for (p in sequence(nsites))
         ylimit = pretty.xylim(u=yrange[,,s],fracexp=0.0,is.log=zlog)
         if (zlog){
            plog     = "y"
            y.at     = pretty.log(ylimit)
            y.labels = sprintf("%g",y.at)
         }else{
            plog     = ""
            y.at     = pretty(ylimit)
            y.labels = sprintf("%g",y.at)
         }#end if
         #---------------------------------------------------------------------------------#


         #------ Set some common features. ------------------------------------------------#
         letitre = paste0(this.desc," - ",sdesc,"\n","Means at equilibrium")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         if (patch.aggr %in% "lorey"){
            lex  = desc.unit(desc="Lorey's height",unit=untab$m)
         }else{
            lex  = desc.unit(desc="Age",unit=untab$yr)
            plog = paste0("x",plog)
         }#end if (patch.aggr %in% "lorey")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ym.patch$default
            fichier = file.path( out.now
                               , paste0("ym_patch-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
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





            #------------------------------------------------------------------------------#
            #      Split the device area into two.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot legend.                                                            #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "bottom"
                  , inset   = 0.0
                  , legend  = paste0(sites$desc," (",toupper(sites$iata),")")
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(nsites)$ncol)
                  , cex     = 0.75
                  , xpd     = TRUE
                  , bty     = "n"
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot simulations.                                                       #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,3.1,1.6))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=plog)
            abline(v=x.at,h=y.at,col=grid.colour,lty="dotted")
            for (p in sequence(nsites)){
               lines(x=xdat[[p]],y=ydat[[p]],lwd=1.5,col=sites$col[p])
            }#end for
            box()
            axis(side=1,at=x.at,labels=x.labels)
            axis(side=2,at=y.at,labels=y.labels,las=1)
            title(main=letitre,xlab=lex,ylab=ley)
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
      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ym.patch)
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
#     Plot the mean seasonal cycle by patch structure.                                     #
#------------------------------------------------------------------------------------------#
if (plot.zm.patch){
   cat0(" + Plot seasonal averages as a function of patch...")


   if (patch.aggr %in% "lorey"){
      x.at     = pretty(x=lorey.at)
      x.labels = sprintf("%g",x.at)
      xlimit   = range(lorey.at)
   }else{
      x.at     = pretty.log(x=age.at)
      x.labels = sprintf("%g",x.at)
      xlimit   = range(age.at)
   }#end if

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.patch        = this.compvar$patchvar
      zlog            = this.compvar$zlog
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't patch-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.patch){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         yrange = array(NA,dim=c(2,nsites,nsimul,nseason.use))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               if (patch.aggr %in% "lorey"){
                  zm.patch   = model[[this.vnam]]$zm.lorey
               }else{
                  zm.patch   = model[[this.vnam]]$zm.age
               }#end if
               if (zlog){
                  zm.patch       = ifelse(zm.patch %>% 0.0,zm.patch,NA)
               }else if (this.vnam %in% "bowen"){
                  zm.patch       = pmax(bmn,pmin(bmx,zm.patch)) + 0. * zm.patch
               }else if (this.vnam %in% "tratio"){
                  zm.patch       = pmax(0.0,pmin(1.0,zm.patch)) + 0. * zm.patch
               }else{
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               for (z in season.use){
                  yrange[,p,s,z] = range(c(yrange[,p,s,z],zm.patch[z,]),finite=TRUE)
               }#end for (z in season.use)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #------ Make the scale for budget variables symmetric. ---------------------------#
         if (this.vnam %in% c("nee","nep","cba")){
            yrange = apply( X      = yrange
                          , MARGIN = c(2,3,4)
                          , FUN    = function(x) c(-1,1)*max(abs(x),na.rm=TRUE)
                          )#end apply
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata           = sites$iata[p]
            this.longname  = sites$desc[p]
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Grab data and arrange them by list.                                     #
            #------------------------------------------------------------------------------#
            xdat = list()
            ydat = list()
            for (s in sequence(nsimul)){
               sname        = simul$name[s]
               model        = res[[iata]][[sname]]
               if (patch.aggr %in% "lorey"){
                  xdat[[s]] = lorey.at
                  zm.patch  = model[[this.vnam]]$zm.lorey
               }else{
                  xdat[[s]] = age.at
                  zm.patch  = model[[this.vnam]]$zm.age
               }#end if
               zm.patch = ifelse(is.finite(zm.patch),zm.patch,NA)
               #---------------------------------------------------------------------------#
               if (this.vnam %in% "bowen"){
                  ydat[[s]] = pmax(bmn,pmin(bmx,zm.patch)) + 0. * zm.patch
               }else if (this.vnam %in% "tratio"){
                  ydat[[s]] = pmax(0.0,pmin(1.0,zm.patch)) + 0. * zm.patch
               }else{
                  ydat[[s]] = zm.patch
               }#end if
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Set some common features. ---------------------------------------------#
            letitre = paste0(this.desc," - ",this.longname)
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            if (patch.aggr %in% "lorey"){
               lex  = desc.unit(desc="Lorey's height",unit=untab$m)
            }else{
               lex  = desc.unit(desc="Age",unit=untab$yr)
               plog = paste0("x",plog)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$zm.patch$variables
               fichier = file.path( out.now
                                  , paste0("zm_patch-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=zsize$width,height=zsize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=zsize$width,height=zsize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=zsize$width,height=zsize$height
                            ,pointsize=ptsz,paper=zsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=zsize$width,height=zsize$height
                     ,pointsize=ptsz,paper=zsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Split the device area into two.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout( mat     =rbind(lo.season$mat.off,rep(1,lo.season$ncol))
                     , heights = c(rep((1.-f.leg)/lo.season$nrow,lo.season$nrow),f.leg)
                     )#end layout
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot legend.                                                         #
               #---------------------------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "center"
                     , inset   = 0.0
                     , legend  = simul$desc
                     , fill    = simul$colour
                     , border  = simul$colour
                     , ncol    = min(3,pretty.box(n=n.sim)$ncol)
                     , title   = expression(bold("Simulation"))
                     , cex     = 0.75
                     , xpd     = TRUE
                     , bty     = "o"
                     )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Plot simulations by season.                                          #
               #---------------------------------------------------------------------------#
               for (z in season.use){
                  ylimit = pretty.xylim(u=yrange[,p,,z],fracexp=0.0,is.log=zlog)
                  if (zlog){
                     plog     = "y"
                     y.at     = pretty.log(ylimit)
                     y.labels = sprintf("%g",y.at)
                  }else{
                     plog     = ""
                     y.at     = pretty(ylimit)
                     y.labels = sprintf("%g",y.at)
                  }#end if
                  if (patch.aggr %in% "age") plog=paste0("x",plog)
                  par(mar=c(4.1,4.6,2.1,0.6))
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                  for (s in sequence(nsimul)){
                     lines(x=xdat[[s]],y=ydat[[s]][z,],lwd=1.5,col=simul$colour[s])
                  }#end for
                  box()
                  axis(side=1,at=x.at,labels=x.labels)
                  axis(side=2,at=y.at,labels=y.labels,las=1)
                  title(main=season.full[z],line=0.8,cex.main=0.8)
               }#end for (z in season.use)
               #---------------------------------------------------------------------------#



               #----- Plot the global title. ----------------------------------------------#
               gtitle( main      = letitre
                     , xlab      = lex
                     , ylab      = ley
                     , off.xlab  = f.leg / (1. + f.leg)
                     , line.main = 2.5
                     , line.ylab = 3.0
                     , cex.axis  = 1.0
                     )#end gtitle
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plot default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname     = simul$name[s]
         sdesc     = simul$desc[s]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Grab data and arrange them by list.                                        #
         #---------------------------------------------------------------------------------#
         xdat = list()
         ydat = list()
         for (p in sequence(nsites)){
            iata         = sites$iata[p]
            longname     = sites$desc[p]
            model        = res[[iata]][[sname]]
            if (patch.aggr %in% "lorey"){
               xdat[[p]] = lorey.at
               zm.patch  = model[[this.vnam]]$zm.lorey
            }else{
               xdat[[p]] = age.at
               zm.patch  = model[[this.vnam]]$zm.age
            }#end if (patch.aggr %in% "lorey")
            zm.patch = ifelse(is.finite(zm.patch),zm.patch,NA)
            if (this.vnam %in% "bowen"){
               ydat[[p]] = pmax(bmn,pmin(bmx,zm.patch)) + 0. * zm.patch
            }else if (this.vnam %in% "tratio"){
               ydat[[p]] = pmax(0.0,pmin(1.0,zm.patch)) + 0. * zm.patch
            }else{
               ydat[[p]] = zm.patch
            }#end if
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #------ Set some common features. ------------------------------------------------#
         letitre = paste0(this.desc," - ",sdesc)
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         if (patch.aggr %in% "lorey"){
            lex  = desc.unit(desc="Lorey's height",unit=untab$m)
         }else{
            lex  = desc.unit(desc="Age",unit=untab$yr)
            plog = paste0("x",plog)
         }#end if (patch.aggr %in% "lorey")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$zm.patch$default
            fichier = file.path( out.now
                               , paste0("zm_patch-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=zsize$width,height=zsize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=zsize$width,height=zsize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=zsize$width,height=zsize$height
                         ,pointsize=ptsz,paper=zsize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=zsize$width,height=zsize$height
                  ,pointsize=ptsz,paper=zsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Split the device area into two.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2.5,0))
            layout( mat     =rbind(lo.season$mat.off,rep(1,lo.season$ncol))
                  , heights = c(rep((1.-f.leg)/lo.season$nrow,lo.season$nrow),f.leg)
                  )#end layout
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot legend.                                                            #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "bottom"
                  , inset   = 0.0
                  , legend  = paste0(sites$desc," (",toupper(sites$iata),")")
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(nsites)$ncol)
                  , title   = expression(bold("Sites"))
                  , cex     = 0.75
                  , xpd     = TRUE
                  , bty     = "o"
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot simulations by season.                                             #
            #------------------------------------------------------------------------------#
            for (z in season.use){
               ylimit = pretty.xylim(u=yrange[,,s,z],fracexp=0.0,is.log=zlog)
               if (zlog){
                  plog     = "y"
                  y.at     = pretty.log(ylimit)
                  y.labels = sprintf("%g",y.at)
               }else{
                  plog     = ""
                  y.at     = pretty(ylimit)
                  y.labels = sprintf("%g",y.at)
               }#end if
               if (patch.aggr %in% "age") plog=paste0("x",plog)
               par(mar=c(4.1,4.6,2.1,0.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=plog)
               for (p in sequence(nsites)){
                  lines(x=xdat[[p]],y=ydat[[p]][z,],lwd=1.5,col=sites$col[p])
               }#end for
               box()
               axis(side=1,at=x.at,labels=x.labels)
               axis(side=2,at=y.at,labels=y.labels,las=1)
               title(main=season.full[z],line=0.8,cex.main=0.8)
            }#end for (z in season.use)
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , xlab      = lex
                  , ylab      = ley
                  , off.xlab  = f.leg / (1. + f.leg)
                  , line.main = 2.5
                  , line.ylab = 3.0
                  , cex.axis  = 1.0
                  )#end gtitle
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
      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ym.patch)
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
#     Plot the mean annual cycle by age structure.                                         #
#------------------------------------------------------------------------------------------#
if (plot.xyz.patch){
   cat0(" + Plot mean annual cycle as a function of patch...")

   xlimit = pretty.xylim(u=c(1.0,13.0),is.log=FALSE)
   if (patch.aggr %in% "lorey"){
      ylog       = FALSE
      y.at       = pretty(x=lorey.at)
      ylimit     = pretty.xylim(lorey.at,is.log=FALSE)
      patch.at   = lorey.at
      n.patch.at = n.lorey.at
   }else{
      ylog       = TRUE
      y.at       = pretty.log(x=age.at)
      ylimit     = pretty.xylim(age.at,is.log=TRUE)
      patch.at   = age.at
      n.patch.at = n.age.at
   }#end if (patch.aggr %in% "lorey")

   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.patch        = this.compvar$patchvar
      zlog            = this.compvar$zlog
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.patch){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         zrange = array(NA,dim=c(2,nsites,nsimul))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               if (patch.aggr %in% "lorey"){
                  mm.patch = model[[this.vnam]]$mm.lorey
               }else{
                  mm.patch = model[[this.vnam]]$mm.age
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               if (zlog){
                  mm.patch = ifelse(mm.patch %>% 0.0,mm.patch,NA)
               }else if (this.vnam %in% "bowen"){
                  mm.patch = pmax(bmn,pmin(bmx,mm.patch)) + 0. * mm.patch
               }else if (this.vnam %in% "tratio"){
                  mm.patch = pmax(0.0,pmin(1.0,mm.patch)) + 0. * mm.patch
               }#end if
               zrange[,p,s] = range(c(zrange[,p,s],mm.patch),finite=TRUE)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #------ Make the scale for budget variables symmetric. ---------------------------#
         if (this.vnam %in% c("nee","nep","cba")){
            zrange = apply( X      = zrange
                          , MARGIN = c(2,3)
                          , FUN    = function(x) c(-1,1)*max(abs(x),na.rm=TRUE)
                          )#end apply
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            zlimit        = pretty.xylim(u=c(zrange[,p,]),is.log=zlog)
            if (zlog){
               z.at       = unique(pretty.log(zlimit,n=ncolours.xyz,forcelog=TRUE))
            }else{
               z.at       = unique(pretty(zlimit,n=ncolours.xyz))
            }#end if
            if (this.vnam %in% c("nee","nep","cba")){
               z.colours  = two.palettes( x     = z.at
                                        , low   = hue.low
                                        , high  = hue.high
                                        , white = 1
                                        , n     = length(z.at)-1
                                        )#end two.palettes
               z.at       = z.colours$breaks
               z.colours  = z.colours$colours
            }else{
               z.colours  = cscheme(n=length(z.at)-1)
            }#end if
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#




            #----- Find the beginning and end of the dry season, and flag them. -----------#
            dryaz    = chron(paste(c(sites$drya[p],sites$dryz[p]),2000,sep="/"))
            xdw      = (numdays(dryaz)-1)/daymax(dryaz) + nummonths(dryaz)
            ydw      = sqrt(prod(ylimit)) + 0 * xdw
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Grab data and arrange them by list.                                     #
            #------------------------------------------------------------------------------#
            xdat = list()
            ydat = list()
            zdat = list()
            x.axis.options = list()
            y.axis.options = list()
            sub.options    = list()
            plot.after     = list()
            for (s in sequence(nsimul)){
               sname               = simul$name[s]
               model               = res[[iata]][[sname]]
               xdat          [[s]] = rep(sequence(12)+0.5,times=n.patch.at)
               ydat          [[s]] = rep(patch.at,each=12)
               if (patch.aggr %in% "lorey"){
                  mm.patch = model[[this.vnam]]$mm.lorey
               }else{
                  mm.patch = model[[this.vnam]]$mm.age
               }#end if
               if (this.vnam %in% "bowen"){
                  zdat       [[s]] = pmax(bmn,pmin(bmx,c(mm.patch)))
               }else if (this.vnam %in% "tratio"){
                  zdat       [[s]] = pmax(0.0,pmin(1.0,c(mm.patch)))
               }else{
                  zdat       [[s]] = c(mm.patch)
               }#end if
               x.axis.options[[s]] = list(side=1,at=1:13
                                         ,labels=substring(c(month.abb,month.abb[1]),1,1))
               y.axis.options[[s]] = list(side=2,las=1,at=y.at)
               sub.options   [[s]] = list(main=simul$desc[s],line=0.6)
               plot.after    [[s]] = list( points = list( x   = xdw
                                                        , y   = ydw
                                                        , pch = 20
                                                        , cex = 0.4
                                                        , col = grey.mg
                                                        )#end list
                                         )#end list
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.desc,this.longname,"Annual means",sep=" - ")
            ley     = desc.unit(desc="Age",unit=untab$yr)
            lex     = "" # "Month"
            lacle   = desc.unit(desc=this.desc,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$xyz.patch$variables
               fichier = file.path( out.now
                                  , paste0("xyz_patch-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=ssize$width,height=ssize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=ssize$width,height=ssize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=ssize$width,height=ssize$height
                            ,pointsize=ptsz,paper=ssize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                     ,pointsize=ptsz,paper=ssize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over simulations.                                               #
               #---------------------------------------------------------------------------#
               image.map ( x               = xdat
                         , y               = ydat
                         , z               = zdat
                         , xlim            = xlimit
                         , ylim            = ylimit
                         , ylog            = ylog
                         , col             = z.colours
                         , levels          = z.at
                         , na.col          = "transparent"
                         , x.axis.options  = x.axis.options
                         , y.axis.options  = y.axis.options
                         , sub.options     = sub.options
                         , main.title      = list( main     = letitre
                                                 , xlab     = lex
                                                 , ylab     = ley
                                                 , cex.main = cex.main
                                                 )#end list
                         , key.title       = list( main     = lacle
                                                 , cex.main = cex.main
                                                 , line     = 1.0
                                                 )#end list
                         , key.log         = zlog
                         , key.vertical    = FALSE
                         , matrix.plot     = TRUE
                         , byrow           = FALSE
                         , plot.after      = plot.after
                         , f.key           = f.leg
                         , smidgen         = 0.04
                         )#end image.map
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plot default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname     = simul$name[s]
         sdesc     = simul$desc[s]
            zlimit        = pretty.xylim(u=c(zrange[,,s]),is.log=zlog)
         if (zlog){
            z.at      = pretty.log(zlimit,n=ncolours.xyz,forcelog=TRUE)
         }else{
            z.at      = pretty(zlimit,n=ncolours.xyz)
         }#end if
         if (this.vnam %in% c("nee","nep","cba")){
            z.colours      = two.palettes( x     = z.at
                                         , low   = hue.low
                                         , high  = hue.high
                                         , white = 1
                                         , n     = length(z.at)-1
                                         )#end two.palettes
            z.at           = z.colours$breaks
            z.colours      = z.colours$colours
         }else{
            z.colours      = cscheme(n=length(z.at)-1)
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Grab data and arrange them by list.                                        #
         #---------------------------------------------------------------------------------#
         xdat = list()
         ydat = list()
         zdat = list()
         x.axis.options = list()
         y.axis.options = list()
         sub.options    = list()
         plot.after     = list()
         for (p in sequence(nsites)){
            iata     = sites$iata[p]
            longname = sites$desc[p]

            #----- Find the beginning and end of the dry season, and flag them. -----------#
            dryaz    = chron(paste(c(sites$drya[p],sites$dryz[p]),2000,sep="/"))
            xdw      = (numdays(dryaz)-1)/daymax(dryaz) + nummonths(dryaz)
            ydw      = sqrt(prod(ylimit)) + 0.* xdw
            #------------------------------------------------------------------------------#

            model               = res[[iata]][[sname]]
            xdat          [[p]] = rep(sequence(12)+0.5,times=n.patch.at)
            ydat          [[p]] = rep(patch.at,each=12)
            if (patch.aggr %in% "lorey"){
               mm.patch         = model[[this.vnam]]$mm.lorey
            }else{
               mm.patch         = model[[this.vnam]]$mm.age
            }#end if (patch.aggr %in% "lorey")


            if (this.vnam %in% "bowen"){
               zdat       [[p]] = pmax(bmn,pmin(bmx,c(mm.patch)))
            }else if (this.vnam %in% "tratio"){
               zdat       [[p]] = pmax(0.0,pmin(1.0,c(mm.patch)))
            }else{
               zdat       [[p]] = c(mm.patch)
            }#end if
            x.axis.options[[p]] = list(side=1,at=1:13
                                      ,labels=substring(c(month.abb,month.abb[1]),1,1))
            y.axis.options[[p]] = list(side=2,las=1,at=y.at)
            sub.options   [[p]] = list(main=longname,line=0.6)
            plot.after    [[p]] = list( abline = list( v   = 1:13
                                                     , h   = y.at
                                                     , col = grid.colour
                                                     , lty = "dotted"
                                                     )#end list
                                      , points = list( x   = xdw
                                                     , y   = ydw
                                                     , pch = 20
                                                     , cex = 0.4
                                                     , col = grey.mg
                                                     )#end list
                                      )#end list
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#


         #------ Set some common features. ------------------------------------------------#
         letitre = paste(this.desc,sdesc,"Annual means",sep=" - ")
         ley     = desc.unit(desc="Age",unit=untab$yr)
         lex     = "" # "Month"
         lacle   = desc.unit(desc=this.desc,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$xyz.patch$default
            fichier = file.path( out.now
                               , paste0("xyz_patch-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=psize$width,height=psize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=psize$width,height=psize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=psize$width*depth,height=psize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height
                         ,pointsize=ptsz,paper=psize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
                  ,pointsize=ptsz,paper=psize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            image.map ( x               = xdat
                      , y               = ydat
                      , z               = zdat
                      , xlim            = xlimit
                      , ylim            = ylimit
                      , ylog            = ylog
                      , col             = z.colours
                      , levels          = z.at
                      , na.col          = "transparent"
                      , x.axis.options  = x.axis.options
                      , y.axis.options  = y.axis.options
                      , sub.options     = sub.options
                      , main.title      = list( main     = letitre
                                              , xlab     = lex
                                              , ylab     = ley
                                              , cex.main = cex.main
                                              )#end list
                      , key.title       = list( main     = lacle
                                              , cex.main = cex.main
                                              , line     = 1.0
                                              )#end list
                      , key.log         = zlog
                      , key.vertical    = FALSE
                      , matrix.plot     = TRUE
                      , byrow           = TRUE
                      , plot.after      = plot.after
                      , f.key           = f.leg
                      , smidgen         = 0.04
                      )#end image.map
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
      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.xyz.patch)
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
#     Plot the mean annual cycle by age structure.                                         #
#------------------------------------------------------------------------------------------#
if (plot.pdf.patch){
   cat0(" + Plot mean annual cycle of variables PDF...")
   xlimit = pretty.xylim(u=c(1.0,13.0),is.log=FALSE)


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high
      is.patch        = this.compvar$patchvar
      #------------------------------------------------------------------------------------#




      #-----  Find the density function. --------------------------------------------------#
      y.levels = seq( from       = min(var.range[1,v,,],na.rm=TRUE)
                    , to         = max(var.range[2,v,,],na.rm=TRUE)
                    , length.out = n.dens
                    )#end seq
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Skip variable if it isn't PFT-dependent.                                       #
      #------------------------------------------------------------------------------------#
      if (is.patch){
         cat0("   - ",this.desc,"...")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         zrange = array(NA,dim=c(2,nsites,nsimul))
         for (p in sequence(nsites)){
            iata     = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               mm.pdf   = model[[this.vnam]]$mm.pdf
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               if (any(is.finite(mm.pdf))){
                  zrange[,p,s] = pdens.range*max(mm.pdf,na.rm=TRUE)
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat0("     * Plotting by sites...")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata           = sites$iata[p]
            this.longname  = sites$desc[p]
            ylimit         = pretty.xylim(u=var.range[,v,p,],fracexp=c(0.05,0.05))
            y.at           = pretty(ylimit)
            y.labels       = sprintf("%g",y.at)
            if (pdens.log){
               z.at           = pretty.log(zrange[,p,],n=ncolours.xyz,forcelog=TRUE)
               z.colours      = cscheme(n=length(z.at)-1)
            }else{
               z.at           = pretty(zrange[,p,],n=ncolours.xyz)
               z.colours      = cscheme(n=length(z.at)-1)
            }#end if
            cat0("       > ",this.longname,"...")
            #------------------------------------------------------------------------------#




            #----- Find the beginning and end of the dry season, and flag them. -----------#
            dryaz    = chron(paste(c(sites$drya[p],sites$dryz[p]),2000,sep="/"))
            xdw      = (numdays(dryaz)-1)/daymax(dryaz) + nummonths(dryaz)
            ydw      = sqrt(prod(ylimit)) + 0 * xdw
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Grab data and arrange them by list.                                     #
            #------------------------------------------------------------------------------#
            xdat = list()
            ydat = list()
            zdat = list()
            x.axis.options = list()
            y.axis.options = list()
            sub.options    = list()
            plot.after     = list()
            for (s in sequence(nsimul)){
               sname               = simul$name[s]
               model               = res[[iata]][[sname]]
               xdat          [[s]] = rep(sequence(12)+0.5,times=n.dens)
               ydat          [[s]] = rep(y.levels,each=12)
               zdat          [[s]] = c(model[[this.vnam]]$mm.pdf)
               x.axis.options[[s]] = list(side=1
                                         ,at=seq(from=1,to=13,by=1)
                                         ,labels=substring(c(month.abb,month.abb[1]),1,1)
                                         ,cex.axis=2.)
               y.axis.options[[s]] = list(side=2,las=1,at=y.at,labels=y.labels,cex.axis=2.)
               sub.options   [[s]] = list(main=simul$desc[s],line=0.6)
               plot.after    [[s]] = list( points = list( x   = xdw
                                                        , y   = ydw
                                                        , pch = 20
                                                        , cex = 0.4
                                                        , col = grey.mg
                                                        )#end list
                                         )#end list
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"(Monthly means)")
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            lex     = "" # "Month"
            lacle   = desc.unit(desc="Probability density function",unit=untab$empty)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$pdf.patch$variables
               fichier = file.path( out.now
                                  , paste0("pdf_patch-",this.vnam,"-",iata,".",outform[o])
                                  )#end file.path
               if (outform[o] %in% "x11"){
                  X11(width=eysize$width,height=eysize$height,pointsize=col.use)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=eysize$width,height=eysize$height,pointsize=col.use)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=eysize$width*depth,height=eysize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=eysize$width*depth,height=eysize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=eysize$width,height=eysize$height
                            ,pointsize=ptsz,paper=eysize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=eysize$width,height=eysize$height
                     ,pointsize=ptsz,paper=eysize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over simulations.                                               #
               #---------------------------------------------------------------------------#
               image.map ( x               = xdat
                         , y               = ydat
                         , z               = zdat
                         , xlim            = xlimit
                         , ylim            = ylimit
                         , col             = z.colours
                         , levels          = z.at
                         , na.col          = "transparent"
                         , x.axis.options  = x.axis.options
                         , y.axis.options  = y.axis.options
                         , sub.options     = sub.options
                         , main.title      = list( main     = ""
                                                 , xlab     = ""
                                                 , ylab     = ley
                                                 , cex.main = cex.main
                                                 )#end list
                         , key.title       = list( main     = lacle
                                                 , cex.main = cex.main
                                                 , line     = 1.0
                                                 )#end list
                         , key.log         = pdens.log
                         , key.vertical    = FALSE
                         , matrix.plot     = FALSE
                         , edge.axes       = TRUE
                         , byrow           = FALSE
                         , plot.after      = plot.after
                         , f.key           = f.leg
                         , smidgen         = 0.04
                         , lo.panel        = pretty.box(n=c(1,n.sim))
                         )#end image.map
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
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         s = sim.default
         cat0("     * Plot default simulation...")
         #----- Get the basic information. ------------------------------------------------#
         sname     = simul$name[s]
         sdesc     = simul$desc[s]
         if (pdens.log){
            z.at      = pretty.log(zrange[,,s],n=ncolours.xyz,forcelog=TRUE)
            z.colours = cscheme(n=length(z.at)-1)
         }else{
            z.at      = pretty(zrange[,,s],n=ncolours.xyz)
            z.colours = cscheme(n=length(z.at)-1)
         }#end if
         ylimit    = pretty.xylim(u=var.range[,v,,s],fracexp=c(0.05,0.05))
         y.at      = pretty(ylimit)
         y.labels  = sprintf("%g",y.at)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Grab data and arrange them by list.                                        #
         #---------------------------------------------------------------------------------#
         xdat = list()
         ydat = list()
         zdat = list()
         x.axis.options = list()
         y.axis.options = list()
         sub.options    = list()
         plot.after     = list()
         for (p in sequence(nsites)){
            iata                = sites$iata[p]
            longname            = sites$desc[p]


            #----- Find the beginning and end of the dry season, and flag them. -----------#
            dryaz    = chron(paste(c(sites$drya[p],sites$dryz[p]),2000,sep="/"))
            xdw      = (numdays(dryaz)-1)/daymax(dryaz) + nummonths(dryaz)
            ydw      = sqrt(prod(ylimit)) + 0.* xdw
            #------------------------------------------------------------------------------#


            model               = res[[iata]][[sname]]
            xdat          [[p]] = rep(sequence(12)+0.5,times=n.dens)
            ydat          [[p]] = rep(y.levels,each=12)
            zdat          [[p]] = c(model[[this.vnam]]$mm.pdf)
            x.axis.options[[p]] = list(side=1,at=1:13
                                      ,labels=substring(c(month.abb,month.abb[1]),1,1))
            y.axis.options[[p]] = list(side=2,las=1,at=y.at,labels=y.labels)
            sub.options   [[p]] = list(main=longname,line=0.6)
            plot.after    [[p]] = list( points = list( x   = xdw
                                                     , y   = ydw
                                                     , pch = 20
                                                     , cex = 0.4
                                                     , col = grey.mg
                                                     )#end list
                                      )#end list
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#


         #------ Set some common features. ------------------------------------------------#
         letitre = paste(sdesc,sep=" - ")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         lex     = "" # "Month"
         lacle   = desc.unit(desc="Probability density function",unit=untab$empty)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$pdf.patch$default
            fichier = file.path( out.now
                               , paste0("pdf_patch-",this.vnam,"-",sname,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=psize$width,height=psize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=psize$width,height=psize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=psize$width*depth,height=psize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height
                         ,pointsize=ptsz,paper=psize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
                  ,pointsize=ptsz,paper=psize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            image.map ( x               = xdat
                      , y               = ydat
                      , z               = zdat
                      , xlim            = xlimit
                      , ylim            = ylimit
                      , col             = z.colours
                      , levels          = z.at
                      , na.col          = "transparent"
                      , x.axis.options  = x.axis.options
                      , y.axis.options  = y.axis.options
                      , sub.options     = sub.options
                      , main.title      = list( main     = letitre
                                              , xlab     = ""
                                              , ylab     = ley
                                              , cex.main = 1.0
                                              )#end list
                      , key.title       = list( main     = lacle
                                              , cex.main = 1.0
                                              , line     = 1.0
                                              )#end list
                      , key.log         = pdens.log
                      , key.vertical    = FALSE
                      , byrow           = TRUE
                      , matrix.plot     = TRUE
                      , plot.after      = plot.after
                      , f.key           = f.leg
                      , smidgen         = 0.04
                      )#end image.map
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
      }#end if (is.szpft)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.pdf.patch)
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
#      Plot the themes.                                                                    #
#------------------------------------------------------------------------------------------#
if (plot.ym.theme){
   cat0(" + Plot annual means for themes...")
   for (th in sequence(nym.theme)){
      this.theme   = ym.theme[[th]]
      theme.vnam   = this.theme$vnam
      theme.desc   = this.theme$desc
      theme.colour = this.theme$colour
      theme.lwd    = this.theme$lwd
      ylog         = this.theme$ylog
      theme.prefix = this.theme$prefix
      theme.title  = this.theme$title
      theme.unit   = this.theme$unit
      cat0("   - ",theme.desc,"...")

      ntheme.vnam  = length(theme.vnam)
      plog         = ifelse(ylog,"y","")


      #----- Get range. -------------------------------------------------------------------#
      xrange = array(data=NA,dim=c(2,nsites,nsimul))
      yrange = array(data=NA,dim=c(2,nsites,nsimul))
      for (p in sequence(nsites)){
         iata = sites$iata[p]
         
         for (s in sequence(nsimul)){
            sname  = simul$name[s]
            model  = res[[iata]][[sname]]
            xrange[,p,s] = range(model$toyear,finite=TRUE)
            for (v in sequence(ntheme.vnam)){
               vnow = model[[theme.vnam[v]]]$ymean
               vnow = ifelse(is.finite(vnow),vnow,NA)
               if (ylog && theme.vnam[v] %in% c("mort","dimort","ncbmort","firemort")){
                  vnow = pmax(0.1,vnow,na.rm=TRUE)
               }#end if
               yrange[,p,s] = range(c(yrange[,p,s],vnow),finite=TRUE)
            }#end for (v in sequence(ntheme.vnam))
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Plot by site (all scenarios).                                                  #
      #------------------------------------------------------------------------------------#
      cat0("     * Plotting by sites...")
      for (p in sequence(nsites)){
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         xlimit        = pretty.xylim(u=xrange[,p,],fracexp=0.0,is.log=FALSE)
         x.at          = pretty(xlimit)
         x.labels      = sprintf("%g",x.at)
         ylimit        = pretty.xylim(u=yrange[,p,],fracexp=0.0,is.log=ylog )
         if (ylog){
            y.at       = pretty.log(ylimit)
            y.labels   = sprintf("%g",y.at)
         }else{
            y.at       = pretty(ylimit)
            y.labels   = sprintf("%g",y.at)
         }#end if
         cat0("       > ",this.longname,"...")
         #---------------------------------------------------------------------------------#




         #------ Set some common features. ------------------------------------------------#
         letitre = paste(this.longname,"Annual means",sep=" - ")
         ley     = desc.unit(desc=theme.title,unit=theme.unit)
         lex     = "Year"
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ym.theme$variables
            fichier = file.path( out.now
                               , paste0("ym_theme-",theme.prefix,"-",iata,".",outform[o])
                               )#end file.path
            if (outform[o] %in% "x11"){
               X11(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if (outform[o] %in% "quartz"){
               quartz(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if(outform[o] %in% "png"){
               png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                  ,pointsize=ptsz,res=depth,bg="transparent")
            }else if(outform[o] %in% "tif"){
               tiff(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                   ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
            }else if(outform[o] %in% "eps"){
               postscript(file=fichier,width=ssize$width,height=ssize$height
                         ,pointsize=ptsz,paper=ssize$paper)
            }else if(outform[o] %in% "pdf"){
               pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                  ,pointsize=ptsz,paper=ssize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2.5,0))
            layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                  , heights = c(rep((1.-f.leg)/lo.simul$nrow,lo.simul$nrow),f.leg)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "center"
                  , inset   = 0.0
                  , legend  = theme.desc
                  , fill    = theme.colour
                  , border  = theme.colour
                  , ncol    = ifelse( ntheme.vnam <= 4
                                    , ntheme.vnam
                                    , min(4,pretty.box(n=ntheme.vnam)$ncol)
                                    )#end ifelse
                  , xpd     = TRUE
                  , bty     = "o"
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               #----- Load variables. -----------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               toyear   = model$toyear
               #---------------------------------------------------------------------------#


               #----- Open window and plot all time series by PFT. ------------------------#
               par(mar=c(4.1,4.6,2.1,1.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=plog)
               abline(h=y.at,v=x.at,col=grid.colour,lty="dotted")
               for (v in rev(sequence(ntheme.vnam))){
                  lines( x    = toyear
                       , y    = model[[theme.vnam[v]]]$ymean
                       , col  = theme.colour[v]
                       , lwd  = theme.lwd   [v]
                       , type = "l"
                       )#end lines
               }#end for (v in rev(ntheme.vnam))
               axis(side=1,las=1,at=x.at,labels=x.labels)
               axis(side=2,las=1,at=y.at,labels=y.labels)
               title(main=simul$desc[s],line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , xlab      = lex
                  , ylab      = ley
                  , line.main = 2.5
                  , off.xlab  = f.leg / (1. + f.leg)
                  )#end gtitle
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
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Plot default simulation.                                                       #
      #------------------------------------------------------------------------------------#
      s = sim.default
      cat0("     * Plotting default simulation...")
      #----- Get the basic information. ---------------------------------------------------#
      sname       = simul$name[s]
      sdesc       = simul$desc[s]
      xlimit      = pretty.xylim(u=xrange[,,s],fracexp=0.0,is.log=FALSE)
      x.at        = pretty(xlimit)
      x.labels    = sprintf("%g",x.at)
      ylimit      = pretty.xylim(u=yrange[,,s],fracexp=0.0,is.log=ylog )
      if (ylog){
         y.at     = pretty.log(ylimit)
         y.labels = sprintf("%g",y.at)
      }else{
         y.at     = pretty(ylimit)
         y.labels = sprintf("%g",y.at)
      }#end if
      #------------------------------------------------------------------------------------#




      #------ Set some common features. ---------------------------------------------------#
      letitre = paste(sdesc,"Annual means",sep=" - ")
      ley     = desc.unit(desc=theme.title,unit=theme.unit)
      lex     = "Year"
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Loop over formats.                                                            #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$ym.theme$default
         fichier = file.path( out.now
                            , paste0("ym_theme-",theme.prefix,"-",sname,".",outform[o])
                            )#end file.path
         if (outform[o] %in% "x11"){
            X11(width=psize$width,height=psize$height,pointsize=col.use)
         }else if (outform[o] %in% "quartz"){
            quartz(width=psize$width,height=psize$height,pointsize=col.use)
         }else if(outform[o] %in% "png"){
            png(filename=fichier,width=psize$width*depth,height=psize$height*depth
               ,pointsize=ptsz,res=depth,bg="transparent")
         }else if(outform[o] %in% "tif"){
            tiff(filename=fichier,width=psize$width*depth,height=psize$height*depth
                ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
         }else if(outform[o] %in% "eps"){
            postscript(file=fichier,width=psize$width,height=psize$height
                      ,pointsize=ptsz,paper=psize$paper)
         }else if(outform[o] %in% "pdf"){
            pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
               ,pointsize=ptsz,paper=psize$paper)
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Split device. -------------------------------------------------------------#
         par(par.user)
         par(oma=c(0,0.25,2.5,0))
         layout( mat     = rbind(lo.site$mat.off,rep(1,times=lo.site$ncol))
               , heights = c(rep((1.-f.leg)/lo.site$nrow,lo.site$nrow),f.leg)
               )#end layout
         #---------------------------------------------------------------------------------#


         #----- Plot legend. --------------------------------------------------------------#
         par(mar=c(0.1,4.6,0.1,2.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x       = "center"
               , inset   = 0.0
               , legend  = theme.desc
               , fill    = theme.colour
               , border  = theme.colour
               , ncol    = ifelse( ntheme.vnam <= 4
                                 , ntheme.vnam
                                 , min(4,pretty.box(n=ntheme.vnam)$ncol)
                                 )#end ifelse
               , xpd     = TRUE
               , bty     = "o"
               )#end legend
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over simulations.                                                     #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            model           = res[[iata]][[sname]]
            toyear          = model$toyear
            #------------------------------------------------------------------------------#


            #----- Open window and plot all time series by PFT. ---------------------------#
            par(mar=c(4.1,4.6,2.1,1.6))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=plog)
            abline(h=y.at,v=x.at,col=grid.colour,lty="dotted")
            for (v in rev(sequence(ntheme.vnam))){
               lines( x    = toyear
                    , y    = model[[theme.vnam[v]]]$ymean
                    , col  = theme.colour[v]
                    , lwd  = theme.lwd   [v]
                    , type = "l"
                    )#end lines
            }#end for (v in rev(ntheme.vnam))
            axis(side=1,las=1,at=x.at,labels=x.labels)
            axis(side=2,las=1,at=y.at,labels=y.labels)
            title(main=this.longname,line=0.5)
            box()
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#



         #----- Plot the global title. ----------------------------------------------------#
         gtitle( main      = letitre
               , xlab      = lex
               , ylab      = ley
               , off.xlab  = f.leg / (1. + f.leg)
               , line.main = 2.5
               , line.ylab = 3.0
               )#end gtitle
         #---------------------------------------------------------------------------------#


         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] %in% c("x11","quartz")){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         dummy = clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
   }#end for (th in sequence(nym.theme))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ym.theme)
#==========================================================================================#
#==========================================================================================#
