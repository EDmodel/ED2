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
srcdir  = "/prj/prjidfca/marcosl/Util/Rsc"      #   Script directory
ibackground    = 0                              # Make figures compatible to background
                                                # 0 -- white
                                                # 1 -- black
                                                # 2 -- dark grey
#----- Output directory -------------------------------------------------------------------#
outroot = file.path(here,paste("longterm_comp_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.hour  = TRUE
reload.range = TRUE
rdata.path   = file.path(here,"RData_longterm")
rdata.suffix = "longterm_ed22.RData"
#------------------------------------------------------------------------------------------#



#----- Default phenology and fire model. --------------------------------------------------#
default.iphen = "phen+02"
default.ifire = "fire03"
emean.yeara   = 1972
emean.yearz   = 2011
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
use.sites  = TRUE
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.struct  = list( name        = c("ble_age30_pft02"  ,"ble_age30_pft05"
                                   ,"sas_age01_pft02"  ,"sas_age01_pft05"
                                   ,"sas_age30_pft02"  ,"sas_age30_pft05"  )
                  , desc        = c("Size 01 + Age 01 + PFT 02"
                                   ,"Size 01 + Age 01 + PFT 05"
                                   ,"Size 80 + Age 01 + PFT 02"
                                   ,"Size 80 + Age 01 + PFT 05"
                                   ,"Size 80 + Age 30 + PFT 02"
                                   ,"Size 80 + Age 30 + PFT 05"
                                   )#end c
                  , verbose     = c("2 PFTs"             ,"5 PFTs"
                                   ,"Size + 2 PFTs"      ,"Size + 5 PFTs"
                                   ,"Size + Age + 2 PFTs","Size + Age + 5 PFTs"
                                   )#end c
                   , colour     = c("#3B24B3","#2996CC"
                                   ,"#990F0F","#E65C17"
                                   ,"#306614","#A3CC52"
                                   )#end c
                   , fgcol      = c("#160959","#0A4766"
                                   ,"#4D0404","#732A06"
                                   ,"#143305","#4B6614"
                                   )#end c
                   , age.interp = c(    FALSE,    FALSE
                                   ,    FALSE,    FALSE
                                   ,     TRUE,     TRUE
                                   )#end c
                   )#end list
#----- List the default simulation. -------------------------------------------------------#
sim.default = 6
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform           = c("pdf")  # Formats for output file.  Supported formats are:
                              #   - "X11" - for printing on screen
                              #   - "eps" - for postscript printing
                              #   - "png" - for PNG printing
                              #   - "pdf" - for PDF printing

byeold            = TRUE      # Remove old files of the given format?

depth             = 96        # PNG resolution, in pixels per inch
paper             = "letter"  # Paper size, to define the plot shape
wpaper            = "legal"   # Wide paper size, to define the plot shape
ptsz              = 22        # Font size.

st.cex.min        = 1.0       # Minimum and maximum sizes for points in the 
st.cex.max        = 2.5       #     Skill and Taylor diagrams
st.lwd.min        = 1.3       # Minimum and maximum sizes for points in the 
st.lwd.max        = 3.0       #     Skill and Taylor diagrams


light.method      = "nls"     # Which method to apply (nls or optim)
light.skew        = FALSE     # Use skew-normal distribution on residuals?
light.n.boot      = 1000      # # of bootstrap iterations (light response fit)
ftnight.n.boot    = 1000      # # of bootstrap iterations (fn mean)

n.quant           = 1024      # # of quantiles to produce the density function.
                              #    We strongly advise to choose a number that is
                              #    a power of two, especially when using EDF 
                              #    (otherwise distributions will be interpolated).
nhour.min         = 16        # Minimum number of hours to use the data.
ust.key.frac      = 1/9       # Fraction of u* plot width used for key
ust.leg.frac      = 1/8       # Fraction of u* plot height used for legend
ncol.ust          = 100       # Number of colours for u*-filter "map"
slz.reference     = -0.50     # Reference level for comparing soil data.
slz.reco.ref      = -0.25     # Reference level for comparing soil data.
n.dens            = 512       # Number of density points. 
dens.min          = 0.00001   # Minimum density (relative to maximum)

keep.gf.low.ustar = FALSE     # Keep data that has been discarded due to low u*, but
                              # otherwise with flux/storage measurements?
                              # TRUE  -- discard only missing data
                              # FALSE -- discard missing data and gap filled due to low
                              #          turbulence.

pft.use           = c(1,2,3,4,16,18)
ncolours.xyz      = 96

#----- Maximum absolute value for Bowen ratio. --------------------------------------------#
bmn = -0.09
bmx =  1.99
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      Colours for model comparison.                                                       #
#------------------------------------------------------------------------------------------#
obs.col   = "#6D6D6D"
obs.fg    = "#131313"
ed22.col  = "#99FF02"
ed22.fg   = "#548901"
alt.col   = "#9D57F8"
alt.fg    = "#2B0071"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Switch controls to plot only the needed ones.                                       #
#------------------------------------------------------------------------------------------#
plot.ts.pft     = c(FALSE,TRUE)[2]
plot.pdf.patch  = c(FALSE,TRUE)[2]
plot.xyz.patch  = c(FALSE,TRUE)[2]
plot.ym.patch   = c(FALSE,TRUE)[2]
plot.ym.theme   = c(FALSE,TRUE)[2]
col.dryseason   = "papayawhip"
col.ust.altern  = "firebrick4"
col.ust.default = "deeppink"
slz.cscheme     = "visible"
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
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = "nplant"
                    , zlog         = TRUE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "ba"
                    , desc         = "Basal area"
                    , unit         = untab$m2om2
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = "nplant"
                    , zlog         = TRUE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "lai"
                    , desc         = "Leaf area index"
                    , unit         = untab$m2lom2
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "gpp"
                    , desc         = "Gross primary productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "reco"
                    , desc         = "Ecosystem respiration"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "nep"
                    , desc         = "Net Ecosystem Productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = FALSE
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
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "wflxca"
                    , desc         = "Water vapour flux"
                    , unit         = untab$kgwom2oday
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "transp"
                    , desc         = "Transpiration"
                    , unit         = untab$kgwom2oday
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = TRUE
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
                    , spider       = FALSE
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
                    , spider       = FALSE
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
                    , spider       = FALSE
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
                    , spider       = FALSE
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
                    , spider       = FALSE
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
                    , spider       = FALSE
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
                    , spider       = FALSE
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
                    , spider       = TRUE
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
                    , spider       = TRUE
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
                    , spider       = TRUE
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
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "can.depth"
                    , desc         = "Mean canopy height"
                    , unit         = untab$m
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "can.area"
                    , desc         = "Mean canopy area"
                    , unit         = untab$empty
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "wood.dens"
                    , desc         = "Mean wood density"
                    , unit         = untab$gocm3
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , szpftvar     = FALSE
                    , patchvar     = TRUE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = FALSE
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
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "ncbmort"
                    , desc         = "Density-dependent mortality"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "dimort"
                    , desc         = "Density independent mortality"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    , spider       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "mort"
                    , desc         = "Mortality rate"
                    , unit         = untab$pcpopoyr
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , szpftvar     = TRUE
                    , patchvar     = FALSE
                    , cohortvar    = FALSE
                    , scalevar     = NA
                    , zlog         = TRUE
                    , spider       = FALSE
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



#----- Find the best set up for plotting all seasons in the same plot. --------------------#
lo.box   = pretty.box(n=nseason-1)
lo.simul = pretty.box(n=nsimul,byrow=FALSE)
lo.site  = pretty.box(n=nsites)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
one.h   = 3.5
one.w   = 3.5 * golden
ssize   = plotsize( proje     = FALSE
                  , stdheight = lo.simul$nrow * 7/6 * one.h
                  , stdwidth  = lo.simul$ncol * one.w
                  )#end plotsize
psize   = plotsize( proje     = FALSE
                  , stdheight = lo.site$nrow  * 7/6 * one.h
                  , stdwidth  = lo.site$ncol  * one.w
                  )#end plotsize
zsize   = plotsize( proje     = FALSE
                  , stdheight = 8.5 * 7/6
                  , stdwidth  = 11
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
   #     Create paths for the time series of annual means.                                 #
   #---------------------------------------------------------------------------------------#
   o.main   = file.path(o.form$main,"ts_pft")
   o.ts.pft = list( main      = o.main
                  , default   = file.path(o.main,"default")
                  , variables = list( main = file.path(o.main,"variables"))
                  )#end list
   if (is.figure){
      if (! file.exists(o.ts.pft$main          )) dir.create(o.ts.pft$main          )
      if (! file.exists(o.ts.pft$default       )) dir.create(o.ts.pft$default       )
      if (! file.exists(o.ts.pft$variables$main)) dir.create(o.ts.pft$variables$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam
      is.szpft         = this.compvar$szpftvar

      #----- Sites. -----------------------------------------------------------------------#
      if (is.szpft){
         o.compvar        = file.path(o.ts.pft$variables$main,this.vnam)
         if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
         o.ts.pft$variables[[this.vnam]] = o.compvar
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$ts.pft = o.ts.pft
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the ym by time of year and patch.                               #
   #---------------------------------------------------------------------------------------#
   o.main  = file.path(o.form$main,"ym_patch")
   o.ym.patch = list( main      = o.main
                 , default   = file.path(o.main,"default")
                 , variables = list( main = file.path(o.main,"variables"))
                 )#end list
   if (is.figure){
      if (! file.exists(o.ym.patch$main          )) dir.create(o.ym.patch$main          )
      if (! file.exists(o.ym.patch$default       )) dir.create(o.ym.patch$default       )
      if (! file.exists(o.ym.patch$variables$main)) dir.create(o.ym.patch$variables$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam
      is.patch         = this.compvar$patchvar

      #----- Sites. -----------------------------------------------------------------------#
      if (is.patch){
         o.compvar        = file.path(o.ym.patch$variables$main,this.vnam)
         if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
         o.ym.patch$variables[[this.vnam]] = o.compvar
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$ym.patch = o.ym.patch
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the PDF by time of year and patch.                               #
   #---------------------------------------------------------------------------------------#
   o.main  = file.path(o.form$main,"pdf_patch")
   o.pdf.patch = list( main      = o.main
                 , default   = file.path(o.main,"default")
                 , variables = list( main = file.path(o.main,"variables"))
                 )#end list
   if (is.figure){
      if (! file.exists(o.pdf.patch$main          )) dir.create(o.pdf.patch$main          )
      if (! file.exists(o.pdf.patch$default       )) dir.create(o.pdf.patch$default       )
      if (! file.exists(o.pdf.patch$variables$main)) dir.create(o.pdf.patch$variables$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam
      is.patch         = this.compvar$patchvar

      #----- Sites. -----------------------------------------------------------------------#
      if (is.patch){
         o.compvar        = file.path(o.pdf.patch$variables$main,this.vnam)
         if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
         o.pdf.patch$variables[[this.vnam]] = o.compvar
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$pdf.patch = o.pdf.patch
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the mean annual cycle by age.                                    #
   #---------------------------------------------------------------------------------------#
   o.main  = file.path(o.form$main,"xyz_patch")
   o.xyz.patch = list( main      = o.main
                 , default   = file.path(o.main,"default")
                 , variables = list( main = file.path(o.main,"variables"))
                 )#end list
   if (is.figure){
      if (! file.exists(o.xyz.patch$main          )) dir.create(o.xyz.patch$main          )
      if (! file.exists(o.xyz.patch$default       )) dir.create(o.xyz.patch$default       )
      if (! file.exists(o.xyz.patch$variables$main)) dir.create(o.xyz.patch$variables$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam
      is.patch         = this.compvar$patchvar

      #----- Sites. -----------------------------------------------------------------------#
      if (is.patch){
         o.compvar        = file.path(o.xyz.patch$variables$main,this.vnam)
         if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
         o.xyz.patch$variables[[this.vnam]] = o.compvar
      }#end if (is.patch)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$xyz.patch = o.xyz.patch
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the theme plots.                                                 #
   #---------------------------------------------------------------------------------------#
   o.main     = file.path(o.form$main,"ym_theme")
   o.ym.theme = list( main      = o.main
                    , default   = file.path(o.main,"default")
                    , variables = list( main = file.path(o.main,"variables"))
                    )#end list
   if (is.figure){
      if (! file.exists(o.ym.theme$main          )) dir.create(o.ym.theme$main          )
      if (! file.exists(o.ym.theme$default       )) dir.create(o.ym.theme$default       )
      if (! file.exists(o.ym.theme$variables$main)) dir.create(o.ym.theme$variables$main)
   }#end if
   for (th in sequence(nym.theme)){
      this.theme       = ym.theme[[th]]
      this.prefix      = this.theme$prefix

      #----- Themes. ----------------------------------------------------------------------#
      o.prefix = file.path(o.ym.theme$variables$main,this.prefix)
      if (is.figure && ! file.exists(o.prefix)) dir.create(o.prefix)
      o.ym.theme$variables[[this.prefix]] = o.prefix
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$ym.theme = o.ym.theme
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the spider web plots.                                            #
   #---------------------------------------------------------------------------------------#
   o.main   = file.path(o.form$main,"spider")
   o.spider = list( main      = o.main
                  , variables = file.path(o.main,"variables")
                  )#end list
   if (is.figure){
      if (! file.exists(o.spider$main     )) dir.create(o.spider$main     )
      if (! file.exists(o.spider$variables)) dir.create(o.spider$variables)
   }#end if
   o.form$spider = o.spider
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
   if (file.exists(rdata.iata)){
      #----- Reload data and copy to the general list. ------------------------------------#
      cat(" + Loading data from ",basename(rdata.iata),"...","\n",sep="")
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
   cat ("   - Loading the ranges...","\n")
   dummy = load(file=rdata.range)
   #---------------------------------------------------------------------------------------#
}else{
   #----- Find ranges. --------------------------------------------------------------------#
   cat ("   - Finding the ranges...","\n")
   age.range = array(data=NA,dim=c(2,nsites,nsimul))
   var.range = array(data=NA,dim=c(2,ncompvar,nsites,nsimul))
   for (p in sequence(nsites)){
      iata          = sites$iata[p]
      longname      = sites$desc[p]
      cat("    > ",sites$desc[p],"...","\n")
      for (s in sequence(nsimul)){
         cat("      # Simulation: ",simul$desc[s],"...","\n")

         #----- Load hourly averages. -----------------------------------------------------#
         ans.name = paste("e",iata,"_wmo_",simul$name[s],"_",default.iphen
                         ,"_",default.ifire,sep="")
         ans.path = file.path(here,ans.name)
         ans.file = file.path(ans.path,"rdata_month",paste(ans.name,".RData",sep=""))
         load(ans.file)
         #---------------------------------------------------------------------------------#



         #----- Grab all relevant times and . ---------------------------------------------#
         esel       = numyears(datum$when) %wr% c(emean.yeara,emean.yearz)
         tomonth    = datum$when[esel]
         nemean     = length(tomonth)
         emean.loop = sequence(nemean)
         agepa      = datum$patch$age
         #---------------------------------------------------------------------------------#



         #----- Update ranges. ------------------------------------------------------------#
         for (e in emean.loop){
            now   = tomonth[e]
            emap  = match(tomonth[e],datum$when)
            mm    = nummonths(now)
            yyyy  = numyears (now)
            stamp = paste("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm),sep="")
            if (simul$age.interp[s]){
               agenow          = pmax(1/12,agepa[[stamp]])
               age.range[,p,s] = range(c(age.range[,p,s],agenow),finite=TRUE)
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
         rm(esel,tomonth,nemean,emean.loop,agepa)
      }#end for (s in sequence(nsimul))
      #------------------------------------------------------------------------------------#
   }#end for (p in loop.sites)
   #---------------------------------------------------------------------------------------#


   #----- Make sure the range is either finite or NA. -------------------------------------#
   var.range = ifelse(is.finite(var.range),var.range,NA)
   dimnames(var.range) = list(c("min","max"),compvar.key,sites$iata,simul$name)
   #---------------------------------------------------------------------------------------#


   #----- Save range. ---------------------------------------------------------------------#
   cat ("   - Finding the ranges...","\n")
   dummy = save(list=c("age.range","var.range"),file=rdata.range)
   #---------------------------------------------------------------------------------------#
}#end if (reload.range && file.exists(rdata.range))
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Define age classes.                                                                  #
#------------------------------------------------------------------------------------------#
age.at   = unique(pretty.log(age.range,n=30,forcelog=TRUE))
n.age.at = length(age.at)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Loop over all sites, variables, and simulations, to prepare the data for the        #
# model comparison.  We want to use only times with actual measurements, so we will        #
# discard model results from times with no observation so all derived quantities have      #
# the same number of defined points (so if measurements are biased towards daytime, the    #
# model will also be equally biased).                                                      #
#------------------------------------------------------------------------------------------#
if (length(loop.sites) != 0) cat (" + Processing missing hourly data...","\n")
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
   cat("   - Site :",this$longname,"...","\n")
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Get all the statistics and actual values for every simulation.                    #
   #---------------------------------------------------------------------------------------#
   cat("    * Aggregate and find statistics for simulations for this site...","\n")
   for (s in sequence(nsimul)){
      cat("      # Simulation: ",simul$desc[s],"...","\n")

      #----- Load hourly averages. --------------------------------------------------------#
      ans.name = paste("e",iata,"_wmo_",simul$name[s],"_",default.iphen,"_",default.ifire
                      ,sep="")
      ans.path = file.path(here,ans.name)
      ans.file = file.path(ans.path,"rdata_month",paste(ans.name,".RData",sep=""))
      load(ans.file)
      #------------------------------------------------------------------------------------#





      #----- Create some variables to describe season and time of the day. ----------------#
      model = list()
      #----- Create time stamp for annual and monthly means. ------------------------------#
      model$toyear  = unique(numyears(datum$when))
      esel          = numyears(datum$when) %wr% c(emean.yeara,emean.yearz)
      model$tomonth = datum$when[esel]
      nymean        = length(model$toyear)
      nemean        = length(model$tomonth)
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #      Load all variables, interpolate them and make the table.                      #
      #------------------------------------------------------------------------------------#
      cat("       ~ Loading variables...","\n")
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
         cat("         > ",this.desc,"...","\n")
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find out how to aggregate the running monthly mean.                        #
         #---------------------------------------------------------------------------------#
         if ( (! is.null(emean)) && (is.null(dim(emean)))){
            emean = matrix(emean,ncol=1)
         }#end if (! is.null(datum$emean))
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Find the annual mean, and crop the monthly mean to the period of interest.  #
         #---------------------------------------------------------------------------------#
         if (! is.null(emean)){
            ymean = qapply( X     = emean
                          , INDEX = numyears(datum$when)
                          , DIM   = 1
                          , FUN   = mean
                          , na.rm = TRUE
                          )#end qapply
            emean = emean[esel,,drop=FALSE]
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Crop the monthly means by PFT and DBH, keeping only the period of interest. #
         #---------------------------------------------------------------------------------#
         if (! is.null(szpft)){
            em.szpft = szpft[esel,,,drop=FALSE]
            ym.szpft = qapply( X     = szpft
                             , INDEX = numyears(datum$when)
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end qapply
         }else{
            em.szpft = NULL
            ym.szpft = NULL
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find mortality rates, "interest" style.                                     #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% c("mort","dimort","ncbmort","firemort")){
            emean    = 100. * ( 1.0 - exp( - emean    ) )
            ymean    = 100. * ( 1.0 - exp( - ymean    ) )
            szpft    = 100. * ( 1.0 - exp( - szpft    ) )
            em.szpft = 100. * ( 1.0 - exp( - em.szpft ) )
            ym.szpft = 100. * ( 1.0 - exp( - ym.szpft ) )
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Aggregate the patch-level data for the last cycle ("equilibrium").          #
         # em.age   -- Equivalent to emean                                                 #
         # mm.age   -- Equivalent to mmean.                                                #
         # pdf.area -- PDF weighted by area and number of occurrences during the period.   #
         #---------------------------------------------------------------------------------#
         if (! is.null(patch)){
            emean.loop = sequence(nemean)
            em.age     = array(data=NA,dim=c(nemean,n.age.at))
            pdf.val    = mapply(FUN=numeric,length=rep(0,times=12),SIMPLIFY=FALSE)
            pdf.wgt    = mapply(FUN=numeric,length=rep(0,times=12),SIMPLIFY=FALSE)
            for (e in emean.loop){
               now     = model$tomonth[e]
               mm      = nummonths(now)
               yyyy    = numyears (now)
               stamp   = paste("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm),sep="")
               vnow    = patch [[stamp]]
               agenow  = pmax(1/12,agepa[[stamp]])
               areanow = areapa[[stamp]]

               #----- Interpolate patch properties to fixed age classes. ------------------#
               if (length(vnow) > 1 && simul$age.interp[s]){
                  age.fun    = splinefun(x=agenow,y=vnow,method="monoH.FC")
                  em.age[e,] = ifelse(age.at %wr% range(agenow),age.fun(u=age.at),NA)
               }else{
                  em.age[e,] = weighted.mean(x=vnow,w=areanow)
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
            mm.age = qapply( X     = em.age
                           , INDEX = nummonths(model$tomonth)
                           , DIM   = 1
                           , FUN   = mean
                           , na.rm = TRUE
                           )#end qapply
            mm.age = ifelse(is.finite(mm.age),mm.age,NA)
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
            em.age = NULL
            mm.age = NULL
            mm.pdf = NULL
         }#end if (! is.null(patch))
         #---------------------------------------------------------------------------------#


         #----- Save variables to a list. -------------------------------------------------#
         model[[this.vnam]] = list( emean    = emean
                                  , ymean    = ymean
                                  , em.szpft = em.szpft
                                  , ym.szpft = ym.szpft
                                  , em.age   = em.age
                                  , mm.age   = mm.age
                                  , mm.pdf   = mm.pdf
                                  )#end list
         rm(list=c("emean","ymean","szpft","em.age","mm.age","mm.pdf"))
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
   cat(" + Saving processed data to ",basename(rdata.iata),"...","\n")
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
#     Plot the time series by PFT for each site.                                           #
#------------------------------------------------------------------------------------------#
if (plot.ts.pft){
   cat(" + Plotting long-term time series...","\n")

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
         cat("   - ",this.desc,"...","\n")


         #---------------------------------------------------------------------------------#
         #     Loop over all sites and simulations, and get the range.                     #
         #---------------------------------------------------------------------------------#
         xrange = array(NA,dim=c(2,nsites,nsimul))
         yrange = array(NA,dim=c(2,nsites,nsimul))
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            for (s in sequence(nsimul)){
               #------ Get the data. ------------------------------------------------------#
               sname    = simul$name[s]
               model    = res[[iata]][[sname]]
               toyear   = model$toyear
               ym.szpft = model[[this.vnam]]$ym.szpft
               nsize    = dim(ym.szpft)[2]
               ym.pft   = ym.szpft[,nsize,]
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               xrange[,p,s] = range(c(xrange[,p,s],toyear)          ,finite=TRUE)
               yrange[,p,s] = range(c(yrange[,p,s],ym.pft[,pft.use]),finite=TRUE)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat("     * Plotting by sites...","\n")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            xlimit          = pretty.xylim(u=xrange[,p,],fracexp=0.0,is.log=FALSE)
            ylimit          = pretty.xylim(u=yrange[,p,],fracexp=0.0,is.log=FALSE)
            cat("       > ",this.longname,"...","\n")
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
               out.now = out[[outform[o]]]$ts.pft$variables[[this.vnam]]
               fichier = file.path( out.now
                                  , paste("ts_pft-",this.vnam,"-",iata,".",outform[o]
                                         ,sep="")
                                  )#end file.path
               if (outform[o] == "x11"){
                  X11(width=ssize$width,height=ssize$height,pointsize=col.use)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=ssize$width,height=ssize$height
                            ,pointsize=ptsz,paper=ssize$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                     ,pointsize=ptsz,paper=ssize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Split device. -------------------------------------------------------#
               par(par.user)
               par(oma=c(0,0,2.5,0))
               layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                     , heights = c(rep(5/lo.simul$nrow,lo.simul$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#


               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,4.6,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "bottom"
                     , inset   = 0.0
                     , legend  = pft$name  [pft.use]
                     , fill    = pft$colour[pft.use]
                     , border  = pft$colour[pft.use]
                     , ncol    = min(3,pretty.box(n=length(pft.use))$ncol)
                     , title   = expression(bold("Plant Functional Type"))
                     , xpd     = TRUE
                     )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Loop over simulations.                                               #
               #---------------------------------------------------------------------------#
               for (s in sequence(nsimul)){
                  #----- Load variables. --------------------------------------------------#
                  sname    = simul$name[s]
                  model    = res[[iata]][[sname]]
                  toyear   = model$toyear
                  ym.szpft = model[[this.vnam]]$ym.szpft
                  nsize    = dim(ym.szpft)[2]
                  ym.pft   = ym.szpft[,nsize,]
                  #------------------------------------------------------------------------#


                  #----- Open window and plot all time series by PFT. ---------------------#
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  grid(col=grid.colour,lty="dotted")
                  for (f in rev(pft.use)){
                     lines( x    = toyear
                          , y    = ym.szpft[,nsize,f]
                          , col  = pft$colour[f]
                          , lwd  = 2.0
                          , type = "l"
                          )#end lines
                  }#end for (f in pft.use)
                  if (lo.simul$bottom[s]) axis(side=1)
                  if (lo.simul$left  [s]) axis(side=2,las=1)
                  title(main=simul$desc[s],line=0.5)
                  box()
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#



               #----- Plot the global title. ----------------------------------------------#
               gtitle( main     = letitre
                     , xlab     = lex
                     , ylab     = ley
                     , off.xlab = 1/12
                     )#end gtitle
               #---------------------------------------------------------------------------#


               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
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
         cat("     * Plotting default simulation...","\n")
         #----- Get the basic information. ------------------------------------------------#
         sname           = simul$name[s]
         sdesc           = simul$desc[s]
         xlimit          = pretty.xylim(u=xrange[,,s],fracexp=0.0,is.log=FALSE)
         ylimit          = pretty.xylim(u=yrange[,,s],fracexp=0.0,is.log=FALSE)
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
                               , paste("ts_pft-",this.vnam,"-",sname,".",outform[o]
                                      ,sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=psize$width,height=psize$height,pointsize=col.use)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height
                         ,pointsize=ptsz,paper=psize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
                  ,pointsize=ptsz,paper=psize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0.25,2.5,0))
            layout( mat     = rbind(lo.site$mat.off,rep(1,times=lo.site$ncol))
                  , heights = c(rep(6/lo.site$nrow,lo.site$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "bottom"
                  , inset   = 0.0
                  , legend  = pft$name  [pft.use]
                  , fill    = pft$colour[pft.use]
                  , border  = pft$colour[pft.use]
                  , ncol    = min(3,pretty.box(n=length(pft.use))$ncol)
                  , title   = expression(bold("Plant Functional Type"))
                  , xpd     = TRUE
                  )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            for (p in sequence(nsites)){
               iata            = sites$iata[p]
               this.longname   = sites$desc[p]
               model           = res[[iata]][[sname]]
               toyear          = model$toyear
               ym.szpft        = model[[this.vnam]]$ym.szpft
               nsize           = dim(ym.szpft)[2]
               ym.pft          = ym.szpft[,nsize,]
               #---------------------------------------------------------------------------#


               #----- Open window and plot all time series by PFT. ------------------------#
               par(mar=lo.site$mar[p,])
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               grid(col=grid.colour,lty="dotted")
               for (f in rev(pft.use)){
                  lines( x    = toyear
                       , y    = ym.szpft[,nsize,f]
                       , col  = pft$colour[f]
                       , lwd  = 2.0
                       , type = "l"
                       )#end lines
               }#end for (f in pft.use)
               if (lo.site$bottom[p]) axis(side=1)
               if (lo.site$left  [p]) axis(side=2,las=1)
               title(main=this.longname,line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main      = letitre
                  , xlab      = lex
                  , ylab      = ley
                  , off.xlab  = 1/12
                  , line.ylab = 3.0
                  )#end gtitle
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
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
#     Plot the mean annual cycle by age structure.                                         #
#------------------------------------------------------------------------------------------#
if (plot.ym.patch){
   cat(" + Plotting equilibrium averages as a function of age...","\n")

   x.at     = pretty.log(x=ceiling(age.at))
   x.labels = sprintf("%g",x.at)
   xlimit   = c( min = max(min(x.at),min(ceiling(age.at)))
               , max = min(max(x.at),max(ceiling(age.at)))
               )#end c

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
         cat("   - ",this.desc,"...","\n")


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
               mm.age   = model[[this.vnam]]$mm.age
               ym.age   = colMeans(mm.age)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               if (zlog){
                  ym.age       = ifelse(ym.age %>% 0.0,ym.age,NA)
                  yrange[,p,s] = range(c(yrange[,p,s],ym.age),finite=TRUE)
               }else if (this.vnam %in% "bowen"){
                  ym.age       = pmax(bmn,pmin(bmx,ym.age))
                  yrange[,p,s] = range(c(yrange[,p,s],ym.age),finite=TRUE)
               }else if (this.vnam %in% "tratio"){
                  ym.age       = pmax(0.0,pmin(1.0,ym.age))
                  yrange[,p,s] = range(c(yrange[,p,s],ym.age),finite=TRUE)
               }else{
                  yrange[,p,s] = range(c(yrange[,p,s],ym.age),finite=TRUE)
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
         cat("     * Plotting by sites...","\n")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata           = sites$iata[p]
            this.longname  = sites$desc[p]
            cat("       > ",this.longname,"...","\n")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Grab data and arrange them by list.                                     #
            #------------------------------------------------------------------------------#
            xdat = list()
            ydat = list()
            for (s in sequence(nsimul)){
               sname        = simul$name[s]
               model        = res[[iata]][[sname]]
               ym.age       = colMeans(model[[this.vnam]]$mm.age)
               xdat   [[s]] = age.at
               if (this.vnam %in% "bowen"){
                  ydat[[s]] = pmax(bmn,pmin(bmx,ym.age))
               }else if (this.vnam %in% "tratio"){
                  ydat[[s]] = pmax(0.0,pmin(1.0,ym.age))
               }else{
                  ydat[[s]] = ym.age
               }#end if
            }#end for (s in sequence(nsimul))
            ylimit = pretty.xylim(u=yrange[,p,],fracexp=0.0,is.log=zlog)
            if (zlog){
               plog     = "xy"
               y.at     = pretty.log(ylimit)
               y.labels = sprintf("%g",y.at)
            }else{
               plog     = "x"
               y.at     = pretty(ylimit)
               y.labels = sprintf("%g",y.at)
            }#end if
            #------------------------------------------------------------------------------#


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.desc," - ",this.longname,"\n","Means at equilibrium"
                           ,sep="")
            lex     = desc.unit(desc="Age",unit=untab$yr)
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ym.patch$variables[[this.vnam]]
               fichier = file.path( out.now
                                  , paste("ym_patch-",this.vnam,"-",iata,".",outform[o]
                                         ,sep="")
                                  )#end file.path
               if (outform[o] == "x11"){
                  X11(width=zsize$width,height=zsize$height,pointsize=col.use)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=zsize$width,height=zsize$height
                            ,pointsize=ptsz,paper=zsize$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=zsize$width,height=zsize$height
                     ,pointsize=ptsz,paper=zsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Split the device area into two.                                      #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(6,1))
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot legend.                                                         #
               #---------------------------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "bottom"
                     , inset   = 0
                     , legend  = simul$desc
                     , fill    = simul$colour
                     , border  = simul$colour
                     , ncol    = min(3,pretty.box(nsimul)$ncol)
                     , title   = expression(bold("Simulations"))
                     , cex     = 0.75
                     , xpd     = TRUE
                     )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Plot simulations.                                                    #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.1,3.1,2.1))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,log=plog)
               abline(v=x.at,h=y.at,col=grid.colour,lty="dotted")
               for (s in sequence(nsimul)){
                  lines(x=xdat[[s]],y=ydat[[s]],lwd=3.0,col=simul$colour[s])
               }#end for
               box()
               axis(side=1,at=x.at,labels=x.labels)
               axis(side=2,at=y.at,labels=y.labels,las=1)
               title(main=letitre,xlab=lex,ylab=ley)
               #---------------------------------------------------------------------------#





               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
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
         cat("     * Plotting default simulation...","\n")
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
            ym.age       = colMeans(model[[this.vnam]]$mm.age)
            xdat   [[p]] = age.at
            if (this.vnam %in% "bowen"){
               ydat[[p]] = pmax(bmn,pmin(bmx,ym.age))
            }else if (this.vnam %in% "tratio"){
               ydat[[p]] = pmax(0.0,pmin(1.0,ym.age))
            }else{
               ydat[[p]] = ym.age
            }#end if
         }#end for (p in sequence(nsites))
         ylimit = pretty.xylim(u=yrange[,,s],fracexp=0.0,is.log=zlog)
         if (zlog){
            plog     = "xy"
            y.at     = pretty.log(ylimit)
            y.labels = sprintf("%g",y.at)
         }else{
            plog     = "x"
            y.at     = pretty(ylimit)
            y.labels = sprintf("%g",y.at)
         }#end if
         #---------------------------------------------------------------------------------#


         #------ Set some common features. ------------------------------------------------#
         letitre = paste(this.desc," - ",sdesc,"\n","Means at equilibrium",sep="")
         lex     = desc.unit(desc="Age",unit=untab$yr)
         lay     = desc.unit(desc=this.desc,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ym.patch$default
            fichier = file.path( out.now
                               , paste("ym_patch-",this.vnam,"-",sname,".",outform[o]
                                      ,sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=zsize$width,height=zsize$height,pointsize=col.use)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=zsize$width,height=zsize$height
                         ,pointsize=ptsz,paper=zsize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=zsize$width,height=zsize$height
                  ,pointsize=ptsz,paper=zsize$paper)
            }#end if
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #      Split the device area into two.                                         #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(6,1))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot legend.                                                            #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "bottom"
                  , inset   = 0
                  , legend  = paste(sites$desc," (",toupper(sites$iata),")",sep="")
                  , fill    = sites$col
                  , border  = sites$col
                  , ncol    = min(3,pretty.box(nsites)$ncol)
                  , title   = expression(bold("Sites"))
                  , cex     = 0.75
                  , xpd     = TRUE
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot simulations.                                                       #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,3.1,2.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,log=plog)
            abline(v=x.at,h=y.at,col=grid.colour,lty="dotted")
            for (p in sequence(nsites)){
               lines(x=xdat[[p]],y=ydat[[p]],lwd=3.0,col=sites$col[p])
            }#end for
            box()
            axis(side=1,at=x.at,labels=x.labels)
            axis(side=2,at=y.at,labels=y.labels,las=1)
            title(main=letitre,xlab=lex,ylab=ley)
            #------------------------------------------------------------------------------#




            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
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
   cat(" + Plotting mean annual cycle as a function of age...","\n")

   y.at   = pretty.log(x=ceiling(age.at))
   xlimit = pretty.xylim(u=c(1.0,13.0),is.log=FALSE)
   ylimit = c( min = max(min(y.at),min(ceiling(age.at)))
             , max = min(max(y.at),max(ceiling(age.at)))
             )#end c

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
         cat("   - ",this.desc,"...","\n")


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
               mm.age   = model[[this.vnam]]$mm.age
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Update range.                                                         #
               #---------------------------------------------------------------------------#
               if (zlog){
                  mm.age       = ifelse(mm.age %>% 0.0,mm.age,NA)
                  zrange[,p,s] = range(c(zrange[,p,s],mm.age),finite=TRUE)
               }else if (this.vnam %in% "bowen"){
                  mm.age       = pmax(bmn,pmin(bmx,mm.age))
                  zrange[,p,s] = range(c(zrange[,p,s],mm.age),finite=TRUE)
               }else if (this.vnam %in% "tratio"){
                  mm.age       = pmax(0.0,pmin(1.0,mm.age))
                  zrange[,p,s] = range(c(zrange[,p,s],mm.age),finite=TRUE)
               }else{
                  zrange[,p,s] = range(c(zrange[,p,s],mm.age),finite=TRUE)
               }#end if
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
         cat("     * Plotting by sites...","\n")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata           = sites$iata[p]
            this.longname  = sites$desc[p]
            if (zlog){
               z.at     = unique(pretty.log(zrange[,p,],n=ncolours.xyz,forcelog=TRUE))
            }else{
               z.at     = unique(pretty(zrange[,p,],n=ncolours.xyz))
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
            cat("       > ",this.longname,"...","\n")
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
               xdat          [[s]] = rep(sequence(12)+0.5,times=n.age.at)
               ydat          [[s]] = rep(age.at,each=12)
               if (this.vnam %in% "bowen"){
                  zdat          [[s]] = pmax(bmn,pmin(bmx,c(model[[this.vnam]]$mm.age)))
               }else if (this.vnam %in% "tratio"){
                  zdat          [[s]] = pmax(0.0,pmin(1.0,c(model[[this.vnam]]$mm.age)))
               }else{
                  zdat          [[s]] = c(model[[this.vnam]]$mm.age)
               }#end if
               x.axis.options[[s]] = list(side=1,at=1:13
                                         ,labels=substring(c(month.abb,month.abb[1]),1,1))
               y.axis.options[[s]] = list(side=2,las=1,at=y.at)
               sub.options   [[s]] = list(main=simul$desc[s],line=0.6)
               plot.after    [[s]] = list( abline = list( v   = 1:13
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
               out.now = out[[outform[o]]]$xyz.patch$variables[[this.vnam]]
               fichier = file.path( out.now
                                  , paste("xyz_patch-",this.vnam,"-",iata,".",outform[o]
                                         ,sep="")
                                  )#end file.path
               if (outform[o] == "x11"){
                  X11(width=ssize$width,height=ssize$height,pointsize=col.use)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=ssize$width,height=ssize$height
                            ,pointsize=ptsz,paper=ssize$paper)
               }else if(outform[o] == "pdf"){
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
                         , ylog            = TRUE
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
                         , f.key           = 1/6
                         , smidgen         = 0.04
                         )#end image.map
               #---------------------------------------------------------------------------#




               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
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
         cat("     * Plotting default simulation...","\n")
         #----- Get the basic information. ------------------------------------------------#
         sname     = simul$name[s]
         sdesc     = simul$desc[s]
         if (zlog){
            z.at      = pretty.log(zrange[,,s],n=ncolours.xyz,forcelog=TRUE)
         }else{
            z.at      = pretty(zrange[,,s],n=ncolours.xyz)
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
            xdat          [[p]] = rep(sequence(12)+0.5,times=n.age.at)
            ydat          [[p]] = rep(age.at,each=12)
            if (this.vnam %in% "bowen"){
               zdat       [[p]] = pmax(bmn,pmin(bmx,c(model[[this.vnam]]$mm.age)))
            }else if (this.vnam %in% "tratio"){
               zdat       [[p]] = pmax(0.0,pmin(1.0,c(model[[this.vnam]]$mm.age)))
            }else{
               zdat       [[p]] = c(model[[this.vnam]]$mm.age)
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
                               , paste("xyz_patch-",this.vnam,"-",sname,".",outform[o]
                                      ,sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=psize$width,height=psize$height,pointsize=col.use)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height
                         ,pointsize=ptsz,paper=psize$paper)
            }else if(outform[o] == "pdf"){
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
                      , ylog            = TRUE
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
                      , f.key           = 1/8
                      , smidgen         = 0.04
                      )#end image.map
            #------------------------------------------------------------------------------#




            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
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
   cat(" + Plotting mean annual cycle of variables PDF...","\n")
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
         cat("   - ",this.desc,"...","\n")


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
                  zrange[,p,s] = c(5.0e-4,1.0)*max(mm.pdf,na.rm=TRUE)
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#
         cat("     * Plotting by sites...","\n")
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata           = sites$iata[p]
            this.longname  = sites$desc[p]
            ylimit         = pretty.xylim(u=var.range[,v,p,],fracexp=c(0.05,0.05))
            y.at           = pretty(ylimit)
            y.labels       = sprintf("%g",y.at)
            z.at           = pretty.log(zrange[,p,],n=ncolours.xyz,forcelog=TRUE)
            z.colours      = cscheme(n=length(z.at)-1)
            cat("       > ",this.longname,"...","\n")
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
                                         ,labels=substring(c(month.abb,month.abb[1]),1,1))
               y.axis.options[[s]] = list(side=2,las=1,at=y.at,labels=y.labels)
               sub.options   [[s]] = list(main=simul$desc[s],line=0.6)
               plot.after    [[s]] = list( abline = list( v   = 1:13
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
            #------------------------------------------------------------------------------#


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.desc,this.longname,"Annual means",sep=" - ")
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            lex     = "" # "Month"
            lacle   = desc.unit(desc="Probability distribution function",unit=untab$empty)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$pdf.patch$variables[[this.vnam]]
               fichier = file.path( out.now
                                  , paste("pdf_patch-",this.vnam,"-",iata,".",outform[o]
                                         ,sep="")
                                  )#end file.path
               if (outform[o] == "x11"){
                  X11(width=ssize$width,height=ssize$height,pointsize=col.use)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=ssize$width,height=ssize$height
                            ,pointsize=ptsz,paper=ssize$paper)
               }else if(outform[o] == "pdf"){
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
                         , col             = z.colours
                         , levels          = z.at
                         , na.col          = "transparent"
                         , x.axis.options  = x.axis.options
                         , y.axis.options  = y.axis.options
                         , sub.options     = sub.options
                         , main.title      = list( main     = letitre
                                                 , xlab     = ""
                                                 , ylab     = ley
                                                 , cex.main = cex.main
                                                 )#end list
                         , key.title       = list( main     = lacle
                                                 , cex.main = cex.main
                                                 , line     = 1.0
                                                 )#end list
                         , key.log         = TRUE
                         , key.vertical    = FALSE
                         , matrix.plot     = TRUE
                         , byrow           = FALSE
                         , plot.after      = plot.after
                         , f.key           = 1/6
                         , smidgen         = 0.04
                         )#end image.map
               #---------------------------------------------------------------------------#




               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
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
         cat("     * Plotting default simulation...","\n")
         #----- Get the basic information. ------------------------------------------------#
         sname     = simul$name[s]
         sdesc     = simul$desc[s]
         z.at      = pretty.log(zrange[,,s],n=ncolours.xyz,forcelog=TRUE)
         z.colours = cscheme(n=length(z.at)-1)
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
         letitre = paste(this.desc,sdesc,sep=" - ")
         ley     = desc.unit(desc=this.desc,unit=this.unit)
         lex     = "" # "Month"
         lacle   = desc.unit(desc="Probability distribution function",unit=untab$empty)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$pdf.patch$default
            fichier = file.path( out.now
                               , paste("pdf_patch-",this.vnam,"-",sname,".",outform[o]
                                      ,sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=psize$width,height=psize$height,pointsize=col.use)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=psize$width*depth,height=psize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=psize$width,height=psize$height
                         ,pointsize=ptsz,paper=psize$paper)
            }else if(outform[o] == "pdf"){
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
                                              , cex.main = cex.main
                                              )#end list
                      , key.title       = list( main     = lacle
                                              , cex.main = cex.main
                                              , line     = 1.0
                                              )#end list
                      , key.log         = TRUE
                      , key.vertical    = FALSE
                      , byrow           = TRUE
                      , matrix.plot     = TRUE
                      , plot.after      = plot.after
                      , f.key           = 1/8
                      , smidgen         = 0.04
                      )#end image.map
            #------------------------------------------------------------------------------#




            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
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
   cat(" + Plotting annual means for themes...","\n")
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
      cat("   - ",theme.desc,"...","\n")

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
      cat("     * Plotting by sites...","\n")
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
         cat("       > ",this.longname,"...","\n")
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
            out.now = out[[outform[o]]]$ym.theme$variables[[theme.prefix]]
            fichier = file.path( out.now
                               , paste("ym_theme-",theme.prefix,"-",iata,".",outform[o]
                                      ,sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=ssize$width,height=ssize$height,pointsize=col.use)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=ssize$width*depth,height=ssize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=ssize$width,height=ssize$height
                         ,pointsize=ptsz,paper=ssize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=ssize$width,height=ssize$height
                  ,pointsize=ptsz,paper=ssize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Split device. ----------------------------------------------------------#
            par(par.user)
            par(oma=c(0,0,2.5,0))
            layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                  , heights = c(rep(5/lo.simul$nrow,lo.simul$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#


            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "bottom"
                  , inset   = 0.0
                  , legend  = theme.desc
                  , fill    = theme.colour
                  , border  = theme.colour
                  , ncol    = ifelse( ntheme.vnam <= 4
                                    , ntheme.vnam
                                    , min(4,pretty.box(n=ntheme.vnam)$ncol)
                                    )#end ifelse
                  , xpd     = TRUE
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
               par(mar=lo.simul$mar[s,])
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
               if (lo.simul$bottom[s]) axis(side=1,at=x.at,labels=x.labels)
               if (lo.simul$left  [s]) axis(side=2,las=1,at=y.at,labels=y.labels)
               title(main=simul$desc[s],line=0.5)
               box()
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#



            #----- Plot the global title. -------------------------------------------------#
            gtitle( main     = letitre
                  , xlab     = lex
                  , ylab     = ley
                  , off.xlab = 1/12
                  )#end gtitle
            #------------------------------------------------------------------------------#


            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
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
      cat("     * Plotting default simulation...","\n")
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
                            , paste("ym_theme-",theme.prefix,"-",sname,".",outform[o]
                                   ,sep="")
                            )#end file.path
         if (outform[o] == "x11"){
            X11(width=psize$width,height=psize$height,pointsize=col.use)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=psize$width*depth,height=psize$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=psize$width,height=psize$height
                      ,pointsize=ptsz,paper=psize$paper)
         }else if(outform[o] == "pdf"){
            pdf(file=fichier,onefile=FALSE,width=psize$width,height=psize$height
               ,pointsize=ptsz,paper=psize$paper)
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Split device. -------------------------------------------------------------#
         par(par.user)
         par(oma=c(0,0.25,2.5,0))
         layout( mat     = rbind(lo.site$mat.off,rep(1,times=lo.site$ncol))
               , heights = c(rep(6/lo.site$nrow,lo.site$nrow),1)
               )#end layout
         #---------------------------------------------------------------------------------#


         #----- Plot legend. --------------------------------------------------------------#
         par(mar=c(0.1,4.6,0.1,2.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x       = "bottom"
               , inset   = 0.0
               , legend  = theme.desc
               , fill    = theme.colour
               , border  = theme.colour
               , ncol    = ifelse( ntheme.vnam <= 4
                                 , ntheme.vnam
                                 , min(4,pretty.box(n=ntheme.vnam)$ncol)
                                 )#end ifelse
               , xpd     = TRUE
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
            par(mar=lo.site$mar[p,])
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
            if (lo.site$bottom[p]) axis(side=1,at=x.at,labels=x.labels)
            if (lo.site$left  [p]) axis(side=2,las=1,at=y.at,labels=y.labels)
            title(main=this.longname,line=0.5)
            box()
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#



         #----- Plot the global title. ----------------------------------------------------#
         gtitle( main      = letitre
               , xlab      = lex
               , ylab      = ley
               , off.xlab  = 1/12
               , line.ylab = 3.0
               )#end gtitle
         #---------------------------------------------------------------------------------#


         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] == "x11"){
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
