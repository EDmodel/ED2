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
outroot = file.path(here,paste("hourly_comp_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.hour  = TRUE
rdata.path   = file.path(here,"RData_often")
rdata.suffix = "hourly_ed22.RData"
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
sites[[n]] = list( iata = "s83"
                 , desc = "Santarem km 83"
                 , pch  =  9
                 , col  = "#E65C17"
                 , fg   = "#732A06"
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
sites[[n]] = list( iata = "cax"
                 , desc = "Caxiuana"
                 , pch  =  0
                 , col  = "#00F3FB"
                 , fg   = "#00AAAF"
                 , drya = "08/19"
                 , dryz = "11/29"
                 )#end list
#use.sites = c("gyf")
#use.sites = c("s67")
#use.sites = c("m34")
#use.sites = c("s83")
#use.sites = c("pdg")
#use.sites = c("s83","pdg","rja","ban")
#use.sites = c("rja","ban")
#use.sites  = c("ban")
use.sites = c("gyf","s67","s83","pdg","rja","m34")
#use.sites = c("gyf","s83","pdg","m34")
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
# sim.struct  = list( name     = c("eft_iust00_ictb00","eft_iust01_ictb00"
#                                 ,"eft_iust00_ictb01","eft_iust01_ictb01"
#                                 ,"eft_iust00_ictb02","eft_iust01_ictb02"
#                                 )#end c
#                   , desc     = c("EFT u*, Leuning"
#                                 ,"Pred. u*, Leuning"
#                                 ,"EFT u*, ED-2.1"
#                                 ,"Pred. u*, ED-2.1"
#                                 ,"EFT u*, Massman"
#                                 ,"Pred. u*, Massman"
#                                 )#end c
#                   , verbose  = c("EFT u*, Leuning"
#                                 ,"Pred. u*, Leuning"
#                                 ,"EFT u*, ED-2.1"
#                                 ,"Pred. u*, ED-2.1"
#                                 ,"EFT u*, Massman"
#                                 ,"Pred. u*, Massman"
#                                 )#end c
#                   , colour   = c("#3B24B3","#2996CC"
#                                 ,"#990F0F","#E65C17"
#                                 ,"#306614","#A3CC52"
#                                 )#end c
#                   , fgcol    = c("#160959","#0A4766"
#                                  ,"#4D0404","#732A06"
#                                  ,"#143305","#4B6614"
#                                  )#end c
#                    )#end list
sim.struct  = list( name     = c("eft_iust00_ictb02","eft_iust01_ictb02")
                  , desc     = c("Prescribed u*"    ,"Predicted u*"     )
                  , verbose  = c("Prescribed u*"    ,"Predicted u*"     )
                  , colour   = c("#3B24B3"          ,"#A3CC52"          )
                  , fgcol    = c("#160959"          ,"#4B6614"          )
                  )#end list
#----- List the default simulation. -------------------------------------------------------#
sim.default = 2
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
ptsz              = 18        # Font size.

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
use.hier.boot     = TRUE      # Use hierarchical sampling?
ust.row.skip      = 1         # When plotting the u-* map, we skip every
                              # ust.row.skip out of ust.row.skip+1 points (to
                              # reduce bad rendering of vector image)
ust.cex.max       = 1.8       # Maximum size for bias.
ust.cex.min       = 0.5       # Maximum size for bias.
bias.max.std      = 3.2       # Maximum bias if the same for all plots
                              # (NULL will find the best for each plot)
n.fit.dmean.min   = 30        # Minimum number of daily averages to fit a curve
n.fit.fmean.min   = 80        # Minimum number of hourly averages to fit a curve

keep.gf.low.ustar = FALSE     # Keep data that has been discarded due to low u*, but
                              # otherwise with flux/storage measurements?
                              # TRUE  -- discard only missing data
                              # FALSE -- discard missing data and gap filled due to low
                              #          turbulence.



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
plot.gpp.light     = c(FALSE,TRUE)[1]
plot.gpp.vpdef     = c(FALSE,TRUE)[1]
plot.gpp.wetness   = c(FALSE,TRUE)[1]
plot.reco.wetness  = c(FALSE,TRUE)[1]
plot.ust.ftnight   = c(FALSE,TRUE)[1]
plot.ust.bias      = c(FALSE,TRUE)[1]
plot.ts.ftnight    = c(FALSE,TRUE)[1]
plot.bp.diel       = c(FALSE,TRUE)[2]
plot.qq.dmean      = c(FALSE,TRUE)[1]
plot.density.dmean = c(FALSE,TRUE)[1]
plot.spider        = c(FALSE,TRUE)[1]
plot.skill.taylor  = c(FALSE,TRUE)[1]
plot.soil.skill    = c(FALSE,TRUE)[1]
make.summ.table    = c(FALSE,TRUE)[1]
density.legend     = FALSE
use.dmean.light    = FALSE
use.dmean.vpdef    = TRUE
use.dmean.wetness  = TRUE
show.dryseason     = TRUE
col.dryseason      = "papayawhip"
col.ust.altern     = "firebrick4"
col.ust.default    = "deeppink"
slz.cscheme        = "visible"
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
n = n + 1
compvar[[ n]] = list( vnam         = "ustar"
                    , symbol       = "u^symbol(\"\\052\")"
                    , desc         = "Friction velocity"
                    , unit         = untab$mos
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#520485"
                    , fg           = "#39025D"
                    , pch          = 9
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "cflxca"
                    , symbol       = "dot(C)[a*e]"
                    , desc         = "Carbon dioxide flux"
                    , unit         = untab$umolcom2os
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , col          = "#46FF32"
                    , fg           = "#31B223"
                    , pch          = 3
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "cflxst"
                    , symbol       = "dot(C)[s*t*o*r]"
                    , desc         = "Carbon dioxide storage"
                    , unit         = untab$umolcom2os
                    , cscheme.mean = "clife"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , col          = "#006715"
                    , fg           = "#00480E"
                    , pch          = 4
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "nee"
                    , symbol       = "N*E*E"
                    , desc         = "Net ecosystem exchange"
                    , unit         = untab$umolcom2os
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , col          = "#46FF32"
                    , fg           = "#31B223"
                    , pch          = 8
                    , ustvar       = TRUE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "nep"
                    , symbol       = "N*E*P"
                    , desc         = "Net ecosystem productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , col          = "#46FF32"
                    , fg           = "#31B223"
                    , pch          = 14
                    , ustvar       = TRUE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "reco"
                    , symbol       = "dot(R)[E*c*o]"
                    , desc         = "Ecosystem respiration"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , col          = "#FF5700"
                    , fg           = "#B23C00"
                    , pch          = 0
                    , ustvar       = TRUE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "gpp"
                    , symbol       = "G*P*P"
                    , desc         = "Gross primary productivity"
                    , unit         = untab$kgcom2oyr
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , col          = "#46FF32"
                    , fg           = "#31B223"
                    , pch          = 2
                    , ustvar       = TRUE
                    , soilvar      = FALSE
                    , sunvar       = TRUE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "parup"
                    , symbol       = "dot(Q)[P*A*R]^symbol(\"\\335\")"
                    , desc         = "Outgoing PAR"
                    , unit         = untab$umolom2os
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#0742C3"
                    , fg           = "#042E88"
                    , pch          = 13
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = TRUE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "rshortup"
                    , symbol       = "dot(Q)[S*W]^symbol(\"\\335\")"
                    , desc         = "Outgoing shortwave radiation"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#FF5700"
                    , fg           = "#B23C00"
                    , pch          = 1
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = TRUE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "rlongup"
                    , symbol       = "dot(Q)[L*W]^symbol(\"\\335\")"
                    , desc         = "Outgoing longwave radiation"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#A00014"
                    , fg           = "#70000E"
                    , pch          = 19
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "hflxca"
                    , symbol       = "dot(Q)[a*e]"
                    , desc         = "Sensible heat flux"
                    , unit         = untab$wom2
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#A00014"
                    , fg           = "#70000E"
                    , pch          = 25
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "wflxca"
                    , symbol       = "dot(W)[a*e]"
                    , desc         = "Water vapour flux"
                    , unit         = untab$kgwom2oday
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , col          = "#0742C3"
                    , fg           = "#042E88"
                    , pch          = 6
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "can.temp"
                    , symbol       = "T[a]"
                    , desc         = "CAS Temperature"
                    , unit         = untab$degC
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#FF5700"
                    , fg           = "#B23C00"
                    , pch          = 15
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "can.shv"
                    , symbol       = "w[a]"
                    , desc         = "CAS specific humidity"
                    , unit         = untab$gwokg
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , col          = "#0742C3"
                    , fg           = "#042E88"
                    , pch          = 12
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam         = "can.co2"
                    , symbol       = "c[a]"
                    , desc         = "CAS CO2 mix. ratio"
                    , unit         = untab$umolcomol
                    , cscheme.mean = "iclife"
                    , hue.low      = "green"
                    , hue.high     = "purple"
                    , col          = "#FF5700"
                    , fg           = "#B23C00"
                    , pch          = 7
                    , ustvar       = FALSE
                    , soilvar      = FALSE
                    , sunvar       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "soil.temp"
                    , symbol       = "T[s*o*i*l]"
                    , desc         = "Soil temperature"
                    , unit         = untab$degC
                    , cscheme.mean = "panoply"
                    , hue.low      = "blue"
                    , hue.high     = "orangered"
                    , col          = "#FF5700"
                    , fg           = "#B23C00"
                    , pch          = 10
                    , ustvar       = FALSE
                    , soilvar      = TRUE
                    , sunvar       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "soil.water"
                    , symbol       = "vartheta[s*o*i*l]"
                    , desc         = "Soil moisture"
                    , unit         = untab$m3wom3
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , col          = "#0742C3"
                    , fg           = "#042E88"
                    , pch          = 5
                    , ustvar       = FALSE
                    , soilvar      = TRUE
                    , sunvar       = FALSE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "soil.wetness"
                    , symbol       = "hat(w)[s*o*i*l]"
                    , desc         = "Soil wetness"
                    , unit         = untab$pc
                    , cscheme.mean = "ipanoply"
                    , hue.low      = "orangered"
                    , hue.high     = "blue"
                    , col          = "#520485"
                    , fg           = "#39025D"
                    , pch          = 18
                    , ustvar       = FALSE
                    , soilvar      = TRUE
                    , sunvar       = FALSE
                    )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Input variables.                                                                     #
#------------------------------------------------------------------------------------------#
n            = 0
control      = list()
n            = n + 1
control[[n]] = list( vnam         = "rshort"
                   , symbol       = "S*W^symbol(\"\\337\")"
                   , desc         = "Incoming shortwave radiation"
                   , unit         = untab$wom2
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#FF5700"
                   , fg           = "#B23C00"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = TRUE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "par"
                   , symbol       = "P*A*R^symbol(\"\\337\")"
                   , desc         = "Incoming PAR"
                   , unit         = untab$umolom2os
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#0742C3"
                   , fg           = "#042E88"
                   , flgvar       = FALSE
                   , soilvar      = FALSE
                   , sunvar       = TRUE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "rlong"
                   , symbol       = "L*W^symbol(\"\\337\")"
                   , desc         = "Incoming longwave radiation"
                   , unit         = untab$wom2
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#A00014"
                   , fg           = "#70000E"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "atm.prss"
                   , symbol       = "p(e*f*t)"
                   , desc         = "Air pressure"
                   , unit         = untab$hpa
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#520485"
                   , fg           = "#39025D"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "atm.temp"
                   , symbol       = "T(e*f*t)"
                   , desc         = "Air temperature"
                   , unit         = untab$degC
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#FF5700"
                   , fg           = "#B23C00"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "atm.shv"
                   , symbol       = "w(e*f*t)"
                   , desc         = "Air specific humidity"
                   , unit         = untab$gwokg
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#0742C3"
                   , fg           = "#042E88"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "atm.vels"
                   , symbol       = "u(e*f*t)"
                   , desc         = "Wind speed"
                   , unit         = untab$mos
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high     = "blue"
                   , col          = "#520485"
                   , fg           = "#39025D"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "rain"
                   , symbol       = "dot(W)(e*f*t)"
                   , desc         = "Precipitation rate"
                   , unit         = untab$kgwom2oday
                   , cscheme.mean = "ipanoply"
                   , hue.low      = "orangered"
                   , hue.high    = "blue"
                   , col          = "#0742C3"
                   , fg           = "#042E88"
                   , flgvar       = TRUE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
n            = n + 1
control[[n]] = list( vnam         = "atm.vpdef"
                   , symbol       = "e(e*f*t)"
                   , desc         = "Vapour pressure deficit"
                   , unit         = untab$hpa
                   , cscheme.mean = "panoply"
                   , hue.low      = "blue"
                   , hue.high     = "orangered"
                   , col          = "#FF5700"
                   , fg           = "#B23C00"
                   , flgvar       = FALSE
                   , soilvar      = FALSE
                   , sunvar       = FALSE
                   )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Statistics.                                                                          #
#------------------------------------------------------------------------------------------#
good = list()
good[[ 1]] = list( vnam      = "bias"
                 , desc      = "Mean bias"
                 , spider    = TRUE
                 , normalise = TRUE
                 )#end list
good[[ 2]] = list( vnam      = "rmse"
                 , desc      = "Root mean square error"
                 , spider    = TRUE
                 , normalise = TRUE
                 )#end list
good[[ 3]] = list( vnam      = "r.squared"
                 , desc      = "Coefficient of determination"
                 , spider    = FALSE
                 , normalise = FALSE
                 )#end list
good[[ 4]] = list( vnam      = "fvue"
                 , desc      = "Fraction of variability unexplained"
                 , spider    = TRUE
                 , normalise = FALSE
                 )#end list
good[[ 5]] = list( vnam      = "sw.stat"
                 , desc      = "Shapiro-Wilk statistic"
                 , spider    = TRUE
                 , normalise = FALSE
                 )#end list
good[[ 6]] = list( vnam      = "ks.stat"
                 , desc      = "Kolmogorov-Smirnov statistic"
                 , spider    = TRUE
                 , normalise = FALSE
                 )#end list
good[[ 7]] = list( vnam      = "lsq.lnlike"
                 , desc      = "Scaled support based on least squares"
                 , spider    = FALSE
                 , normalise = FALSE
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
simcol.key     = sim.struct$colour
simfgc.key     = sim.struct$fgcol
simlty.key     = rep("solid",times=n.sim)
simcex.key     = rep(2.0    ,times=n.sim)
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
control.key  = list.2.data.frame(control)$vnam
compvar.key  = list.2.data.frame(compvar)$vnam
compvar.pch  = list.2.data.frame(compvar)$pch
compvar.sym  = parse(text=list.2.data.frame(compvar)$symbol)
good.key     = list.2.data.frame(good   )$vnam
season.key   = season.list
diel.key     = c("night"  ,"rise.set","day"    ,"all.hrs"  ,"dmean"     
                ,"fnmean")
diel.desc    = c("Night"  ,"Rise/Set","Day"    ,"All hours","Daily mean"
                ,"Semi-monthly mean")
diel.default = c(TRUE     ,FALSE     ,TRUE     ,TRUE       ,FALSE       
                ,FALSE    )
diel.col     = c("#520485","#FF5700" ,"#A00014"  ,"#46FF32","#0742C3"   
                ,"#F5C858")
moment.key   = c("mean","variance","skewness","kurtosis")
moment.desc  = c("Mean","Variance","Skewness","Kurtosis")
hour.num     = seq(from=0,to=23,by=3)
hour.key     = c( paste( sprintf("%2.2i",(hour.num - 1) %% 24)
                      , sprintf("%2.2i",(hour.num + 2) %% 24)
                      , sep = "."
                      )#end paste
               )#end c
hour.label  = paste( sprintf("%2.2i",(hour.num - 1) %% 24)
                   , sprintf("%2.2i",(hour.num + 2) %% 24)
                   , sep = "-"
                   )#end paste
hour.desc   = paste(hour.label,"UTC",sep = " ")
#------------------------------------------------------------------------------------------#



#----- Set the various dimensions associated with variables, simulations, and sites. ------#
nsites     = length(sites$iata )
nsimul     = length(simul.key  )
ncompvar   = length(compvar.key)
ncontrol   = length(control.key)
ngood      = length(good.key   )
nseason    = length(season.key )
ndiel      = length(diel.key   )
nhour      = length(hour.key   )
nmoment    = length(moment.key )
#------------------------------------------------------------------------------------------#





#----- Save special diel times. -----------------------------------------------------------#
diel.all.hrs = max(c(0,match("all.hrs",diel.key)),na.rm=TRUE)
diel.dmean   = max(c(0,match("dmean"  ,diel.key)),na.rm=TRUE)
diel.fnmean  = max(c(0,match("fnmean" ,diel.key)),na.rm=TRUE)
#------------------------------------------------------------------------------------------#



season.suffix = paste(sprintf("%2.2i",sequence(nseason)),tolower(season.key),sep="-")


#----- Load observations. -----------------------------------------------------------------#
obser.file = paste(srcdir,"LBA_MIP.nogapfill.RData",sep="/")
load(file=obser.file)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size    = plotsize(proje=FALSE,paper=paper )
wsize   = plotsize(proje=FALSE,paper=wpaper)
tsize   = plotsize(proje=FALSE,stdheight=wsize$height*8/7,stdwidth=wsize$width*7/6)
#------------------------------------------------------------------------------------------#



#----- Find the best set up for plotting all seasons in the same plot. --------------------#
lo.box   = pretty.box(n=nseason-1)
lo.simul = pretty.box(n=nsimul,byrow=FALSE)
lo.site  = pretty.box(n=nsites)
#------------------------------------------------------------------------------------------#






#------ Set some common features for fortnightly means. -----------------------------------#
fnmean.year   = 2004
fnmean.when   = fnyear.2.chron(fortnight = sequence(yr.ftnight),year = fnmean.year)
fnmean.axis   = pretty.ftnight(fnmean.when)
fnmean.at     = fnmean.axis$levels
fnmean.limit  = range(fnmean.at)
fnmean.labels = fnmean.axis$labels
dfnmean.when  = 0.2 * mean(diff(fnmean.when))
#------------------------------------------------------------------------------------------#


#------ Find the minimum number of hours to consider some statistics. ---------------------#
ef.min = 1.0 # nhour.min / day.hr
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
   #     Create paths for the "spider web" plots.                                          #
   #---------------------------------------------------------------------------------------#
   o.spider              = list()
   o.spider$main         = file.path(o.form$main  ,"spider"           )
   o.spider$default.var  = file.path(o.spider$main,"default_variables")
   o.spider$default.site = file.path(o.spider$main,"default_sites"    )
   if (is.figure && ! file.exists(o.spider$main        )) dir.create(o.spider$main        )
   if (is.figure && ! file.exists(o.spider$default.var )) dir.create(o.spider$default.var )
   if (is.figure && ! file.exists(o.spider$default.site)) dir.create(o.spider$default.site)
   for (d in sequence(ndiel)){
      this.diel        = diel.key [d]
      o.diel           = list()
      o.diel$main      = file.path(o.spider$main,this.diel  )
      o.diel$sites     = file.path(o.diel$main  ,"sites"    )
      o.diel$variables = file.path(o.diel$main  ,"variables")
      if (is.figure){
         if (! file.exists(o.diel$main     )) dir.create(o.diel$main     )
         if (! file.exists(o.diel$sites    )) dir.create(o.diel$sites    )
         if (! file.exists(o.diel$variables)) dir.create(o.diel$variables)
      }#end if (is.figure)
      #------------------------------------------------------------------------------------#
      o.spider[[this.diel]] = o.diel
   }#end for (d in 1:ndiel)
   o.form$spider = o.spider
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of fortnightly means.                            #
   #---------------------------------------------------------------------------------------#
   o.main = file.path(o.form$main,"ts_ftnight")
   o.ts.ftnight = list( main    = o.main
                      , default = file.path(o.main,"default")
                      , sites   = list( main = file.path(o.main,"sites"     )  )
                      )#end list
   if (is.figure){
      if (! file.exists(o.ts.ftnight$main      )) dir.create(o.ts.ftnight$main      )
      if (! file.exists(o.ts.ftnight$default   )) dir.create(o.ts.ftnight$default   )
      if (! file.exists(o.ts.ftnight$sites$main)) dir.create(o.ts.ftnight$sites$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam

      #----- Sites. -----------------------------------------------------------------------#
      o.compvar        = file.path(o.ts.ftnight$sites$main,this.vnam)
      if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
      o.ts.ftnight$sites[[this.vnam]] = o.compvar
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$ts.ftnight = o.ts.ftnight
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of variables that are function of u* filter.     #
   #---------------------------------------------------------------------------------------#
   o.main        = file.path(o.form$main,"ust_ftnight")
   o.ust.ftnight = list( main    = o.main
                       , default = file.path(o.main,"default")
                       , sites   = list( main = file.path(o.main,"sites"     )  )
                       )#end list
   if (is.figure){
      if (! file.exists(o.ust.ftnight$main      )) dir.create(o.ust.ftnight$main      )
      if (! file.exists(o.ust.ftnight$default   )) dir.create(o.ust.ftnight$default   )
      if (! file.exists(o.ust.ftnight$sites$main)) dir.create(o.ust.ftnight$sites$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam
      this.ust         = this.compvar$vnam %in% c("nee","nep","gpp","reco")

      #----- Sites. -----------------------------------------------------------------------#
      o.compvar        = file.path(o.ust.ftnight$sites$main,this.vnam)
      if (is.figure && this.ust && ! file.exists(o.compvar)) dir.create(o.compvar)
      o.ust.ftnight$sites[[this.vnam]] = o.compvar
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$ust.ftnight = o.ust.ftnight
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the time series of variables that are function of u* filter.     #
   #---------------------------------------------------------------------------------------#
   o.main     = file.path(o.form$main,"ust_bias")
   o.ust.bias = list( main    = o.main
                    , default = file.path(o.main,"default")
                    , sites   = list( main = file.path(o.main,"sites"     )  )
                    )#end list
   if (is.figure){
      if (! file.exists(o.ust.bias$main      )) dir.create(o.ust.bias$main      )
      if (! file.exists(o.ust.bias$default   )) dir.create(o.ust.bias$default   )
      if (! file.exists(o.ust.bias$sites$main)) dir.create(o.ust.bias$sites$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam
      this.ust         = this.compvar$vnam %in% c("nee","nep","gpp","reco")

      #----- Sites. -----------------------------------------------------------------------#
      o.compvar        = file.path(o.ust.bias$sites$main,this.vnam)
      if (is.figure && this.ust && ! file.exists(o.compvar)) dir.create(o.compvar)
      o.ust.bias$sites[[this.vnam]] = o.compvar
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$ust.bias = o.ust.bias
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the box plots of the mean diel.                                  #
   #---------------------------------------------------------------------------------------#
   o.main = file.path(o.form$main,"bp_diel")
   o.bp.diel = list( main    = o.main
                   , default = file.path(o.main,"default")
                   , sites   = list( main = file.path(o.main,"sites"))
                   )#end list
   if (is.figure){
      if (! file.exists(o.bp.diel$main      )) dir.create(o.bp.diel$main      )
      if (! file.exists(o.bp.diel$default   )) dir.create(o.bp.diel$default   )
      if (! file.exists(o.bp.diel$sites$main)) dir.create(o.bp.diel$sites$main)
   }#end if
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam

      #----- Sites. -----------------------------------------------------------------------#
      o.compvar        = file.path(o.bp.diel$sites$main,this.vnam)
      if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
      o.bp.diel$sites[[this.vnam]] = o.compvar
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   o.form$bp.diel = o.bp.diel
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the Q-Q plots of daily means.                                    #
   #---------------------------------------------------------------------------------------#
   o.qq.dmean         = list()
   o.qq.dmean$main    = file.path(o.form$main    ,"qq_dmean")
   o.qq.dmean$sites   = file.path(o.qq.dmean$main,"sites"   )
   o.qq.dmean$default = file.path(o.qq.dmean$main,"default" )
   if (is.figure && ! file.exists(o.qq.dmean$main)){
      dir.create(o.qq.dmean$main)
   }#end if
   if (is.figure && ! file.exists(o.qq.dmean$sites  )){
      dir.create(o.qq.dmean$sites)
   }#end if
   if (is.figure && ! file.exists(o.qq.dmean$default)){
      dir.create(o.qq.dmean$default)
   }#end if
   o.form$qq.dmean  = o.qq.dmean
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the density function of daily means.                             #
   #---------------------------------------------------------------------------------------#
   o.density.dmean          = list()
   o.density.dmean$main     = file.path(o.form$main         ,"density_dmean")
   o.density.dmean$sites    = file.path(o.density.dmean$main,"sites"        )
   o.density.dmean$default  = file.path(o.density.dmean$main,"default"      )
   if (is.figure && ! file.exists(o.density.dmean$main)){
      dir.create(o.density.dmean$main)
   }#end if
   if (is.figure && ! file.exists(o.density.dmean$sites)){
      dir.create(o.density.dmean$sites)
   }#end if
   if (is.figure && ! file.exists(o.density.dmean$default)){
      dir.create(o.density.dmean$default)
   }#end if
   o.form$density.dmean  = o.density.dmean
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the skill diagrams.                                              #
   #---------------------------------------------------------------------------------------#
   o.skill              = list( main = file.path(o.form$main,"skill"))
   o.skill$default.var  = file.path(o.skill$main ,"default_variables")
   o.skill$default.site = file.path(o.skill$main ,"default_sites")
   if (is.figure && ! file.exists(o.skill$main        )) dir.create(o.skill$main        )
   if (is.figure && ! file.exists(o.skill$default.var )) dir.create(o.skill$default.var )
   if (is.figure && ! file.exists(o.skill$default.site)) dir.create(o.skill$default.site)
   for (d in sequence(ndiel)){
      this.diel            = diel.key [d]
      o.diel               = file.path(o.skill$main,this.diel)
      if (is.figure && ! file.exists(o.diel)) dir.create(o.diel)
      o.skill[[this.diel]] = o.diel
      #------------------------------------------------------------------------------------#
   }#end for (d in sequence(ndiel))
   o.form$skill = o.skill
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the Taylor diagrams.                                             #
   #---------------------------------------------------------------------------------------#
   o.taylor              = list( main = file.path(o.form$main,"taylor"))
   o.taylor$default.var  = file.path(o.taylor$main ,"default_variables")
   o.taylor$default.site = file.path(o.taylor$main ,"default_sites")
   if (is.figure && ! file.exists(o.taylor$main        )) dir.create(o.taylor$main        )
   if (is.figure && ! file.exists(o.taylor$default.var )) dir.create(o.taylor$default.var )
   if (is.figure && ! file.exists(o.taylor$default.site)) dir.create(o.taylor$default.site)
   for (d in sequence(ndiel)){
      this.diel        = diel.key [d]
      o.diel           = file.path(o.taylor$main,this.diel)
      if (is.figure && ! file.exists(o.diel)) dir.create(o.diel)
      o.taylor[[this.diel]] = o.diel
      #------------------------------------------------------------------------------------#
   }#end for (d in sequence(ndiel))
   o.form$taylor = o.taylor
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the skill diagrams.                                              #
   #---------------------------------------------------------------------------------------#
   o.soil.skill           = list( main = file.path(o.form$main,"soil_skill"))
   o.soil.skill$variables = file.path(o.soil.skill$main,"variables")
   o.soil.skill$default   = file.path(o.soil.skill$main,"default"  )
   if(is.figure && ! file.exists(o.soil.skill$main     )) dir.create(o.soil.skill$main     )
   if(is.figure && ! file.exists(o.soil.skill$variables)) dir.create(o.soil.skill$variables)
   if(is.figure && ! file.exists(o.soil.skill$default  )) dir.create(o.soil.skill$default  )
   o.form$soil.skill = o.soil.skill
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the light response curve.                                        #
   #---------------------------------------------------------------------------------------#
   o.light         = list()
   o.light$main    = file.path(o.form$main ,"light"  )
   o.light$sites   = file.path(o.light$main,"sites"  )
   o.light$default = file.path(o.light$main,"default")
   if (is.figure && ! file.exists(o.light$main   )) dir.create(o.light$main   )
   if (is.figure && ! file.exists(o.light$sites  )) dir.create(o.light$sites  )
   if (is.figure && ! file.exists(o.light$default)) dir.create(o.light$default)
   o.form$light  = o.light
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the VPD response curve.                                          #
   #---------------------------------------------------------------------------------------#
   o.vpdef         = list()
   o.vpdef$main    = file.path(o.form$main ,"vpdef")
   o.vpdef$sites   = file.path(o.vpdef$main,"sites"  )
   o.vpdef$default = file.path(o.vpdef$main,"default")
   if (is.figure && ! file.exists(o.vpdef$main   )) dir.create(o.vpdef$main   )
   if (is.figure && ! file.exists(o.vpdef$sites  )) dir.create(o.vpdef$sites  )
   if (is.figure && ! file.exists(o.vpdef$default)) dir.create(o.vpdef$default)
   o.form$vpdef  = o.vpdef
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the soil wetness response curve.                                 #
   #---------------------------------------------------------------------------------------#
   o.wetness         = list()
   o.wetness$main    = file.path(o.form$main   ,"wetness")
   o.wetness$sites   = file.path(o.wetness$main,"sites"  )
   o.wetness$default = file.path(o.wetness$main,"default")
   if (is.figure && ! file.exists(o.wetness$main   )) dir.create(o.wetness$main   )
   if (is.figure && ! file.exists(o.wetness$sites  )) dir.create(o.wetness$sites  )
   if (is.figure && ! file.exists(o.wetness$default)) dir.create(o.wetness$default)
   o.form$wetness  = o.wetness
   #---------------------------------------------------------------------------------------#


   #----- Save the full list to the main path list. ---------------------------------------#
   out[[this.form]] = o.form
   #---------------------------------------------------------------------------------------#
}#end for (o in 1:nout)
#------------------------------------------------------------------------------------------#







#------------------------------------------------------------------------------------------#
#      Retrieve all data that already exists.                                              #
#------------------------------------------------------------------------------------------#
eft        = list()
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
      eft.iata    = paste("eft",iata,sep=".")
      res.iata    = paste("res",iata,sep=".")
      eft[[iata]] = get(eft.iata)
      res[[iata]] = get(res.iata)
      rm(eft.iata,res.iata)
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
#      Loop over all sites, variables, and simulations, to prepare the data for the        #
# model comparison.  We want to use only times with actual measurements, so we will        #
# discard model results from times with no observation so all derived quantities have      #
# the same number of defined points (so if measurements are biased towards daytime, the    #
# model will also be equally biased).                                                      #
#------------------------------------------------------------------------------------------#
if (length(loop.sites) == 0) cat (" + Processing missing hourly data...","\n")
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
   this$sim      = list()
   this$ans      = list()
   cat("   - Site :",this$longname,"...","\n")
   #---------------------------------------------------------------------------------------#


   #----- Grab the observation. -----------------------------------------------------------#
   obser      = get(paste("obs",iata,sep="."))
   #---------------------------------------------------------------------------------------#



   #----- Make some dimension info. -------------------------------------------------------#
   obser$nwhen     = length(obser$when      )
   obser$nustar    = length(obser$ust.filter)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     We are about to go through the simulations.  We now save the time span and soil   #
   # depths from the original simulation, we will trim the observations to the simulation  #
   # period and depth.                                                                     #
   #---------------------------------------------------------------------------------------#
   nobser     = obser$nwhen
   slza.obser = min(obser$slz)
   #---------------------------------------------------------------------------------------#



   #----- Load hourly averages. -----------------------------------------------------------#
   ans.name   = paste("t",iata,"_",simul$name[1],sep="")
   ans.path   = file.path(here,ans.name)
   ans.file   = file.path(ans.path,"rdata_hour",paste(ans.name,".RData",sep=""))
   load(ans.file)
   nmodel     = length(model$when)
   slza.model = min(model$slz)
   #---------------------------------------------------------------------------------------#




   #=======================================================================================#
   #=======================================================================================#
   #      In case the models and observations don't match perfectly, we trim observations. #
   # This is done only once, otherwise the models don't have the same length, which is     #
   # unacceptable.                                                                         #
   #---------------------------------------------------------------------------------------#
   if (nobser != nmodel || slza.obser < slza.model){
      cat("     * Resolving mismatches between model and observations dimensions...","\n")

      #----- Get the model time range. ----------------------------------------------------#
      this.whena = min(as.numeric(model$when))
      this.whenz = max(as.numeric(model$when))
      sel.when   = obser$when >= this.whena & obser$when <= this.whenz
      #------------------------------------------------------------------------------------#



      #----- Crop depths that are deeper, but not shallower. ------------------------------#
      sel.slz   = obser$slz >= slza.model
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over all observed variables and trim them to the same model period.       #
      #------------------------------------------------------------------------------------#
      for (vname in names(obser)){
         #----- Inquire variable dimensions. ----------------------------------------------#
         dim.obs       = dim   (obser[[vname]])
         len.obs       = length(obser[[vname]])
         is.this.mat   = ! is.null(dim.obs)
         is.this.chron = is.chron(obser[[vname]])
         if (is.this.mat){
            is.this.when = dim.obs[1] == nobser
            is.this.slz  = dim.obs[2] == obser$nzg
         }else{
            is.this.when = len.obs    == nobser
            is.this.slz  = len.obs    == obser$nzg
         }#end if (is.mat)
         #---------------------------------------------------------------------------------#


         #----- Crop times. ---------------------------------------------------------------#
         if (is.this.mat && is.this.when){
            obser[[vname]] = obser[[vname]][sel.when,]
         }else if (is.this.when){
            obser[[vname]] = obser[[vname]][sel.when ]
         }else if (is.this.chron){
            sel.now = obser[[vname]] >= this.whena & obser[[vname]] <= this.whenz
            obser[[vname]] = obser[[vname]][sel.now]
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Crop soil levels. ---------------------------------------------------------#
         if (is.this.slz && is.this.mat){
            obser[[vname]] = obser[[vname]][,sel.slz]
         }else if (is.this.slz){
            obser[[vname]] = obser[[vname]][ sel.slz]
         }#end if (is.this.slz && is.this.mat)
         #---------------------------------------------------------------------------------#
      }#end for (vname in names(obser))
      #------------------------------------------------------------------------------------#
   }#end if (nobser != nmodel || slza.obser < slza.model)
   rm (model)
   #---------------------------------------------------------------------------------------#



   #------ Since times and depths may have changed, update dimensions. --------------------#
   obser$nwhen = length(obser$when)
   obser$nzg   = length(obser$slz )
   nobser      = obser$nwhen
   slza.obser  = min(obser$slz)
   #---------------------------------------------------------------------------------------#



   #----- Create some variables to describe season and time of the day. -------------------#
   obser$season    = season(obser$when,add.year=FALSE)
   obser$diel      = 1 + (! obser$nighttime) + obser$highsun
   obser$today     = dates(obser$when)
   obser$fortnight = numfortnight  (obser$when)
   obser$toftnight = fnyear.2.chron(obser$when)
   obser$hr.idx    = 3 * floor(obser$hour / 3)
   #---------------------------------------------------------------------------------------#




   #------  Define unique dates and seasons for daily and fortnightly means. --------------#
   obser$dmean.when    = unique(obser$today)
   obser$fnmean.when   = fnmean.when
   obser$dmean.season  = season(when = obser$dmean.when , add.year=FALSE)
   obser$fnmean.season = season(when = obser$fnmean.when, add.year=FALSE)
   obser$ndmean        = length(obser$dmean.when )
   obser$nfnmean       = length(obser$fnmean.when)
   #---------------------------------------------------------------------------------------#



   #----- Make some labels and suffixes for soil variables. -------------------------------#
   obser$slz.key     = paste("_slz",sprintf("%+06.2f",obser$slz),sep="")
   obser$slz.desc    = paste(" (z=",sprintf("%6.2f",abs(obser$slz)),"m)",sep="")
   #---------------------------------------------------------------------------------------#



   #----- Make some two-dimensional variables for soil and u* plots. ----------------------#
   obser$soil.fnmean.when = matrix( data  = obser$fnmean.when
                                  , nrow  = obser$nfnmean
                                  , ncol  = obser$nzg
                                  , byrow = FALSE
                                  )#end matrix
   obser$ust.fnmean.when  = matrix( data  = obser$fnmean.when
                                  , nrow  = obser$nfnmean
                                  , ncol  = obser$nustar
                                  , byrow = FALSE
                                  )#end matrix
   obser$soil.slz         = matrix( data  = obser$slz
                                  , nrow  = obser$nfnmean
                                  , ncol  = obser$nzg
                                  , byrow = TRUE
                                  )#end matrix
   obser$ust.ustmin       = matrix( data  = obser$ust.filter
                                  , nrow  = obser$nfnmean
                                  , ncol  = obser$nustar
                                  , byrow = TRUE
                                  )#end matrix
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Fix the measured flags for GPP and ecosystem respiration.                        #
   #---------------------------------------------------------------------------------------#
   cat("     * Fix flags for ecosystem estimates...","\n")
   #----- Flags to tell that CO2 flux and storage were measured. --------------------------#
   if (all(! obser$measured.cflxst)){
      measured.co2         = obser$measured.cflxca
   }else{
      measured.co2         = obser$measured.cflxca & obser$measured.cflxst
   }#end if
   ust.measured.co2        = matrix( data     = measured.co2
                                   , nrow     = obser$nwhen
                                   , ncol     = obser$nustar
                                   , dimnames = dimnames(obser$ust.nee)
                                   )#end matrix
   #----- Decide what to keep depending on whether to keep low turbulence periods. --------#
   if (keep.gf.low.ustar){
      #----- Keep all measured times, even if it was filled due to low turbulence. --------#
      measured.eco     = measured.co2
      alt.measured.eco = measured.co2
      ust.measured.eco = ust.measured.co2
      #------------------------------------------------------------------------------------#
   }else{
      #----- Keep only measurements with enough turbulence. -------------------------------#
      measured.eco         = obser$measured.nee
      alt.measured.eco     = obser$alt.measured.nee
      #------------------------------------------------------------------------------------#


      #----- Create a quick filter for different u*. --------------------------------------#
      ustar.eco        = ifelse(obser$measured.ustar & measured.co2,obser$ustar,NA)
      ust.measured.eco = outer(ustar.eco,obser$ust.filter,FUN='%>=%')
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #------ Replace the flags of the ecosystem variables. ----------------------------------#
   obser$measured.nep      = measured.eco
   obser$measured.nee      = measured.eco
   obser$measured.gpp      = measured.eco
   obser$measured.reco     = measured.eco
   obser$alt.measured.nep  = alt.measured.eco
   obser$alt.measured.nee  = alt.measured.eco
   obser$alt.measured.gpp  = alt.measured.eco
   obser$alt.measured.reco = alt.measured.eco
   obser$ust.measured.nep  = ust.measured.eco
   obser$ust.measured.nee  = ust.measured.eco
   obser$ust.measured.gpp  = ust.measured.eco
   obser$ust.measured.reco = ust.measured.eco
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      Fix the measured flags for radiation.                                            #
   #---------------------------------------------------------------------------------------#
   cat("     * Fix flags for radiation variables...","\n")
   #----- Assumed that solar radiation is measured if any component is measured. ----------#
   measured.sun          = obser$measured.par | obser$measured.rshort
   obser$measured.par    = measured.sun
   obser$measured.rshort = measured.sun



   #----- In case soil.temp is soil.tmp, replace the name. --------------------------------#
   if ("soil.tmp" %in% names(obser)){
      ist  = match("soil.tmp",names(obser))
      imst = match("measured.soil.tmp",names(obser))
      names(obser)[c(ist,imst)] = c("soil.temp","measured.soil.temp")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- In case can.temp is soil.tmp, replace the name. ---------------------------------#
   if ("can.tmp" %in% names(obser)){
      ist  = match("can.tmp",names(obser))
      imst = match("measured.can.tmp",names(obser))
      names(obser)[c(ist,imst)] = c("can.temp","measured.can.temp")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Compute soil wetness for each soil layer. ---------------------------------------#
   obser$soil.wetness          = apply(X=obser$soil.water,MARGIN=2,FUN=percentil,trim=0.02)
   obser$measured.soil.wetness = obser$measured.soil.water
   #---------------------------------------------------------------------------------------#



   #----- Compute vapour pressure deficit. ------------------------------------------------#
   obser$atm.vpdef          = 0.01 * vpdefil( pres = obser$atm.prss * 100.
                                            , temp = obser$atm.temp + t00
                                            , humi = obser$atm.shv  * 0.001
                                            )#end vpdefil
   obser$measured.atm.vpdef = pmax(obser$measured.atm.shv,obser$measured.atm.temp)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Longwave radiation is lost if the model resumes run, and it is not the only one   #
   # lost.  Throw away all observations when it happens.                                   #
   #---------------------------------------------------------------------------------------#
   discard = rep(FALSE,times=obser$nwhen)
   for (s in sequence(nsimul)){
      ans.name   = paste("t",iata,"_",simul$name[s],sep="")
      ans.path   = file.path(here,ans.name)
      ans.file   = file.path(ans.path,"rdata_hour",paste(ans.name,".RData",sep=""))
      load(ans.file)
      discard = discard | model$rlongup < 1.0
      rm(model)
   }#end for (s in sequence(nsimul))
   for (v in sequence(ncompvar+ncontrol)){
      #----- Load information. ------------------------------------------------------------#
      if (v <= ncompvar){
         this.compvar  = compvar[[v]]
      }else{
         this.compvar  = control[[v-ncompvar]]
      }#end if
      this.vnam     = this.compvar$vnam
      this.measured = paste("measured",this.vnam,sep=".")
      this.soilvar  = this.compvar$soilvar

      if (this.soilvar){
         obser[[this.measured]][discard,] = FALSE
      }else{
         obser[[this.measured]][discard ] = FALSE
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Count the number of variables that have been measured, and use a backward filter  #
   #---------------------------------------------------------------------------------------#
   tot.measured = rep(0,times=obser$nwhen)
   nflags                 = 0
   for (v in sequence(ncontrol)){
      this.compvar  = control[[v]]
      this.measured = paste("measured",this.vnam,sep=".")
      this.vnam     = this.compvar$vnam
      this.flgvar   = this.compvar$flgvar
      if (this.flgvar){
         nflags        = nflags + 1
         tot.measured  = tot.measured + obser[[this.measured]]
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   obser$driver.score = filter( x        = tot.measured
                              , filter   = rep(1/day.hr,times=day.hr)
                              , method   = "convolution"
                              , sides    = 1L
                              , circular = TRUE
                              )#end filter
   obser$driver.use   = obser$driver.score >= max(0.,nflags - 1.5)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the daily mean and create flags on whether to use it or not.                #
   #---------------------------------------------------------------------------------------#
   cat("    * Aggregating observations...","\n")
   for (v in sequence(ncompvar+ncontrol)){
      #----- Load information. ------------------------------------------------------------#
      if (v <= ncompvar){
         this.compvar  = compvar[[v]]
      }else{
         this.compvar  = control[[v-ncompvar]]
      }#end if
      this.vnam     = this.compvar$vnam
      this.dmean    = paste("dmean"   ,this.vnam,sep=".")
      this.fnmean   = paste("fnmean"  ,this.vnam,sep=".")
      this.fnq025   = paste("fnq025"  ,this.vnam,sep=".")
      this.fnq975   = paste("fnq975"  ,this.vnam,sep=".")
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = this.compvar$desc
      this.unit     = this.compvar$unit
      this.soilvar  = this.compvar$soilvar
      this.sunvar   = this.compvar$sunvar
      cat("      # ",this.desc,"...","\n")
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Check whether this is a soil variable.                                         #
      #------------------------------------------------------------------------------------#
      if (this.soilvar){
         driver.ok = matrix( data = rep(obser$driver.use,times=obser$nzg)
                           , nrow = obser$nwhen
                           , ncol = obser$nzg
                           )#end matrix
      }else{
         driver.ok = obser$driver.use
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Discard all gap-filled entries. ----------------------------------------------#
      keep = is.finite(obser[[this.vnam]]) & obser[[this.measured]] & driver.ok
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Replace gap-filled values by NA, and find the daily mean.  We replace the      #
      # original values by these, so we can eliminate model values easily.                 #
      #------------------------------------------------------------------------------------#
      obs.now     = ifelse(keep,obser[[this.vnam]],NA)
      if (! is.matrix(obs.now)){
         obs.now = matrix(obs.now,ncol=1,dimnames=list(names(obs.now),"42.0"))
      }#end if (! is.matrix(obs.now))
      dmean.obser = qapply( X     = obs.now
                          , INDEX = obser$today
                          , DIM   = 1
                          , FUN   = mean
                          , na.rm = TRUE
                          )#end tapply
      dmean.count = qapply( X     = is.finite(obs.now)
                          , INDEX = obser$today
                          , DIM   = 1
                          , FUN   = sum
                          , na.rm = TRUE
                          )#end tapply
      dmean.obser = ifelse(dmean.count >= nhour.min,dmean.obser,NA)
      #------------------------------------------------------------------------------------#



      #----- Save number of layers. -------------------------------------------------------#
      nlyrs = ncol(obs.now)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Solve fortnightly means and find confidence intervals for the mean.            #
      #------------------------------------------------------------------------------------#
      #----- For loop over columns (sorry, Steve...) --------------------------------------#
      fnmean.obser = NULL
      fnq025.obser = NULL
      fnq975.obser = NULL
      for (cc in sequence(nlyrs)){

         keep  = is.finite(obs.now[,cc])

         if (any(keep)){
            #----- Find mean fortnightly period. ------------------------------------------#
            data.in    = data.frame( x         = obs.now[,cc]
                                   , fortnight = obser$fortnight
                                   , year      = obser$year
                                   , hour      = obser$hr.idx
                                   )#end data.frame
            if (use.hier.boot){
               fn.boot = bhier.fortnight.mean( data.in = data.in
                                             , R       = ftnight.n.boot
                                             , ci      = 0.95
                                             )#end bhier.fortnight.mean
            }else{
               fn.boot = boot.fortnight.mean ( data.in = data.in
                                             , R       = ftnight.n.boot
                                             , ci      = 0.95
                                             )#end boot.fortnight.mean
            }#end if
            fnmean.now = fn.boot$expected
            fnq025.now = fn.boot$qlow
            fnq975.now = fn.boot$qhigh
            #------------------------------------------------------------------------------#
         }else{
            #----- Empty data. ------------------------------------------------------------#
            fnmean.now = rep(x=NA,times=obser$nfnmean)
            fnq025.now = rep(x=NA,times=obser$nfnmean)
            fnq975.now = rep(x=NA,times=obser$nfnmean)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Save this column. ---------------------------------------------------------#
         fnmean.obser = cbind(fnmean.obser,fnmean.now)
         fnq025.obser = cbind(fnq025.obser,fnq025.now)
         fnq975.obser = cbind(fnq975.obser,fnq975.now)
         rm(fnmean.now,fnq025.now,fnq975.now)
         #---------------------------------------------------------------------------------#
      }#end for (cc in sequence(ncol(obs.now)))
      #------------------------------------------------------------------------------------#




      #----- Save the data. ---------------------------------------------------------------#
      obser[[this.vnam  ]] = obs.now
      obser[[this.dmean ]] = dmean.obser
      obser[[this.fnmean]] = fnmean.obser
      obser[[this.fnq025]] = fnq025.obser
      obser[[this.fnq975]] = fnq975.obser
      #------------------------------------------------------------------------------------#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#






      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #     Check whether to find alternative variables.                                   #
      #------------------------------------------------------------------------------------#
      this.alt      = paste("alt"           ,this.vnam,sep=".")
      this.dmean    = paste("alt","dmean"   ,this.vnam,sep=".")
      this.fnmean   = paste("alt","fnmean"  ,this.vnam,sep=".")
      this.fnq025   = paste("alt","fnq025"  ,this.vnam,sep=".")
      this.fnq975   = paste("alt","fnq975"  ,this.vnam,sep=".")
      this.measured = paste("alt","measured",this.vnam,sep=".")
      if (this.alt %in% names(obser)){
         cat("        > Averages for alternative data set...","\n")
         #----- Discard all gap-filled entries. -------------------------------------------#
         keep = is.finite(obser[[this.alt]]) & obser[[this.measured]] & driver.ok
         #---------------------------------------------------------------------------------#

         obs.now     = ifelse(keep,obser[[this.alt]],NA)
         if (! is.matrix(obs.now)){
            obs.now = matrix(obs.now,ncol=1,dimnames=list(names(obs.now),"42.0"))
         }#end if (! is.matrix(obs.now))
         dmean.obser = qapply( X     = obs.now
                             , INDEX = obser$today
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end tapply
         dmean.count = qapply( X     = is.finite(obs.now)
                             , INDEX = obser$today
                             , DIM   = 1
                             , FUN   = sum
                             , na.rm = TRUE
                             )#end tapply
         dmean.obser = ifelse(dmean.count >= nhour.min,dmean.obser,NA)
         #---------------------------------------------------------------------------------#



         #----- Save number of layers. ----------------------------------------------------#
         nlyrs = ncol(obs.now)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Solve fortnightly means and find confidence intervals for the mean.         #
         #---------------------------------------------------------------------------------#
         #----- For loop over columns (sorry, Steve...) -----------------------------------#
         fnmean.obser = NULL
         fnq025.obser = NULL
         fnq975.obser = NULL
         for (cc in sequence(nlyrs)){

            keep  = is.finite(obs.now[,cc])

            if (any(keep)){
               #----- Find mean fortnightly period. ---------------------------------------#
               data.in    = data.frame( x         = obs.now[,cc]
                                      , fortnight = obser$fortnight
                                      , year      = obser$year
                                      , hour      = obser$hr.idx
                                      )#end data.frame
               if (use.hier.boot){
                  fn.boot = bhier.fortnight.mean( data.in = data.in
                                                , R       = ftnight.n.boot
                                                , ci      = 0.95
                                                )#end bhier.fortnight.mean
               }else{
                  fn.boot = boot.fortnight.mean ( data.in = data.in
                                                , R       = ftnight.n.boot
                                                , ci      = 0.95
                                                )#end boot.fortnight.mean
               }#end if
               fnmean.now = fn.boot$expected
               fnq025.now = fn.boot$qlow
               fnq975.now = fn.boot$qhigh
               #---------------------------------------------------------------------------#
            }else{
               #----- Empty data. ---------------------------------------------------------#
               fnmean.now = rep(x=NA,times=obser$nfnmean)
               fnq025.now = rep(x=NA,times=obser$nfnmean)
               fnq975.now = rep(x=NA,times=obser$nfnmean)
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#


            #----- Save this column. ------------------------------------------------------#
            fnmean.obser = cbind(fnmean.obser,fnmean.now)
            fnq025.obser = cbind(fnq025.obser,fnq025.now)
            fnq975.obser = cbind(fnq975.obser,fnq975.now)
            rm(fnmean.now,fnq025.now,fnq975.now)
            #------------------------------------------------------------------------------#
         }#end for (cc in sequence(ncol(obs.now)))
         #---------------------------------------------------------------------------------#




         #----- Save the data. ------------------------------------------------------------#
         obser[[this.alt   ]] = obs.now
         obser[[this.dmean ]] = dmean.obser
         obser[[this.fnmean]] = fnmean.obser
         obser[[this.fnq025]] = fnq025.obser
         obser[[this.fnq975]] = fnq975.obser
         #---------------------------------------------------------------------------------#
      }#end if
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#





      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #     Check whether to find u*-filter dependent variables.                           #
      #------------------------------------------------------------------------------------#
      this.ust      = paste("ust"           ,this.vnam,sep=".")
      this.dmean    = paste("ust","dmean"   ,this.vnam,sep=".")
      this.fnmean   = paste("ust","fnmean"  ,this.vnam,sep=".")
      this.measured = paste("ust","measured",this.vnam,sep=".")
      if (this.ust %in% names(obser)){
         cat("        > Averages for u*-matrix...","\n")


         #----- Make matrix for driver check. ---------------------------------------------#
         driver.ok = matrix( data = rep(obser$driver.use,times=obser$nustar)
                           , nrow = obser$nwhen
                           , ncol = obser$nustar
                           )#end matrix
         #---------------------------------------------------------------------------------#


         #----- Discard all gap-filled entries. -------------------------------------------#
         keep = is.finite(obser[[this.ust]]) & obser[[this.measured]] & driver.ok
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Replace gap-filled values by NA, and find the daily mean.  We replace the   #
         # original values by these, so we can eliminate model values easily.              #
         #---------------------------------------------------------------------------------#
         obs.now     = ifelse(keep,obser[[this.ust]],NA)
         dmean.obser = qapply( X     = obs.now
                             , INDEX = obser$today
                             , DIM   = 1
                             , FUN   = mean
                             , na.rm = TRUE
                             )#end tapply
         dmean.count = qapply( X     = is.finite(obs.now)
                             , INDEX = obser$today
                             , DIM   = 1
                             , FUN   = sum
                             , na.rm = TRUE
                             )#end tapply
         dmean.obser = ifelse(dmean.count >= nhour.min,dmean.obser,NA)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Solve fortnightly means.                                                    #
         #---------------------------------------------------------------------------------#
         if (sum(is.finite(c(obs.now))) > 0){

            #----- Save number of years. --------------------------------------------------#
            year.num = sort(unique(obser$year))
            n.year   = length(year.num)
            #------------------------------------------------------------------------------#



            #----- Find mean fortnightly period. ------------------------------------------#
            fnmean.1st = qapply( X      = obs.now
                               , INDEX  = list(obser$hr.idx,obser$fortnight,obser$year)
                               , DIM    = 1
                               , FUN    = mean
                               , na.rm  = TRUE
                               )#end apply
            hh = match(as.numeric(dimnames(fnmean.1st)[[1]]),hour.num            )
            ff = match(as.numeric(dimnames(fnmean.1st)[[2]]),sequence(yr.ftnight))
            yy = match(as.numeric(dimnames(fnmean.1st)[[3]]),year.num            )
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Collapse fortnightly periods.  Make sure all periods are                 #
            # defined in the structure, if none of them were selected, make them           #
            #NA.                                                                           #
            #------------------------------------------------------------------------------#
            fnmean.obser  = array ( data     = NA
                                  , dim      = c(nhour,yr.ftnight,n.year,obser$nustar)
                                  , dimnames = list( hour.num
                                                   , sequence(yr.ftnight)
                                                   , year.num
                                                   , obser$ust.filter
                                                   )#end list
                                  )#end array
            fnmean.obser[hh,ff,yy,] = fnmean.1st
            fnmean.obser = apply(X=fnmean.obser,MARGIN=c(2,3,4),FUN=mean)
            fnmean.obser = apply(X=fnmean.obser,MARGIN=c(1,3)  ,FUN=mean,na.rm=TRUE)
            if (this.vnam %in% c("gpp","reco")){
               fnmean.obser = ifelse( is.finite(fnmean.obser) & fnmean.obser > 0.
                                    , fnmean.obser
                                    , NA
                                    )#end ifelse
            }else{
               fnmean.obser = ifelse( is.finite(fnmean.obser),fnmean.obser,NA)
            }#end if
            #------------------------------------------------------------------------------#

         }else{
            #----- Leave everything empty. ------------------------------------------------#
            fnmean.obser = matrix( data     = NA
                                 , nrow     = yr.ftnight
                                 , ncol     = obser$nustar
                                 , dimnames = list(sequence(yr.ftnight),obser$ust.filter)
                                 )#end fnmean.obser
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Save the data. ------------------------------------------------------------#
         obser[[this.ust   ]] = obs.now
         obser[[this.dmean ]] = dmean.obser
         obser[[this.fnmean]] = fnmean.obser
         #---------------------------------------------------------------------------------#
      }#end if (this.ust %in% names(obser))
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Get all the statistics and actual values for every simulation.                    #
   #---------------------------------------------------------------------------------------#
   cat("    * Aggregate and find statistics for simulations for this site...","\n")
   for (s in sequence(nsimul)){
      cat("      # Simulation: ",simul$desc[s],"...","\n")

      #----- Load hourly averages. --------------------------------------------------------#
      ans.name   = paste("t",iata,"_",simul$name[s],sep="")
      ans.path   = file.path(here,ans.name)
      ans.file   = file.path(ans.path,"rdata_hour",paste(ans.name,".RData",sep=""))
      load(ans.file)
      nmodel     = length(model$when)
      slza.model = min(model$slz)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Make sure the model and the observations cover the same period and soil depth, #
      # and that times are synchronised.                                                   #
      #------------------------------------------------------------------------------------#
      if (nobser != nmodel || slza.obser < slza.model){
         cat(" -> Simulation:"   ,ans.name          ,"\n")
         cat(" -> Length(obser):",length(obser$when),"\n")
         cat(" -> Length(model):",length(model$when),"\n")
         cat(" -> SLZA(obser):"  ,slza.obser        ,"\n")
         cat(" -> SLZA(model):"  ,slza.model        ,"\n")
         stop(" Model and obser must have the same length")
      }else if (any(as.numeric(model$when-obser$when) > 1/48,na.rm=TRUE)){
         stop(" All times in the model and observations must match!!!")
      }#end if
      #------------------------------------------------------------------------------------#





      #----- Create some variables to describe season and time of the day. ----------------#
      model$nwhen         = obser$nwhen
      model$season        = obser$season
      model$diel          = obser$diel
      model$today         = obser$today
      model$fortnight     = obser$fortnight
      model$toftnight     = obser$toftnight
      model$hr.idx        = obser$hr.idx
      model$dmean.when    = obser$dmean.when
      model$fnmean.when   = obser$fnmean.when
      model$dmean.season  = obser$dmean.season
      model$fnmean.season = obser$fnmean.season
      model$ndmean        = obser$ndmean
      model$nfnmean       = obser$nfnmean
      model$year          = obser$year
      model$mon           = obser$mon
      model$day           = obser$day
      model$hour          = obser$hour
      model$nustar        = obser$nustar
      model$ust.filter    = obser$ust.filter
      #------------------------------------------------------------------------------------#





      #----- In case atm.temp is atm.tmp, replace the name. -------------------------------#
      if ("atm.tmp" %in% names(model)){
         ist  = match("atm.tmp",names(model))
         names(model)[ist] = c("atm.temp")
      }#end if
      #------------------------------------------------------------------------------------#





      #----- In case soil.temp is soil.tmp, replace the name. -----------------------------#
      if ("soil.tmp" %in% names(model)){
         ist  = match("soil.tmp",names(model))
         names(model)[ist] = c("soil.temp")
      }#end if
      #------------------------------------------------------------------------------------#





      #----- In case can.temp is can.tmp, replace the name. -------------------------------#
      if ("can.tmp" %in% names(model)){
         ist  = match("can.tmp",names(model))
         names(model)[ist] = c("can.temp")
      }#end if
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #      Interpolate soil variables to the levels with data.                           #
      #------------------------------------------------------------------------------------#
      cat("       ~ Interpolating soil variables...","\n")
      for (v in sequence(ncompvar+ncontrol)){
         #----- Load information. ---------------------------------------------------------#
         if (v <= ncompvar){
            this.compvar  = compvar[[v]]
         }else{
            this.compvar  = control[[v-ncompvar]]
         }#end if
         this.vnam     = this.compvar$vnam
         this.dmean    = paste("dmean"   ,this.vnam,sep=".")
         this.fnmean   = paste("fnmean"  ,this.vnam,sep=".")
         this.fnq025   = paste("fnq025"  ,this.vnam,sep=".")
         this.fnq975   = paste("fnq975"  ,this.vnam,sep=".")
         this.measured = paste("measured",this.vnam,sep=".")
         this.desc     = this.compvar$desc
         this.unit     = this.compvar$unit
         this.soilvar  = this.compvar$soilvar
         this.sunvar   = this.compvar$sunvar
         #---------------------------------------------------------------------------------#


         #----- Make sure this is really a soil variable. ---------------------------------#
         this.soilvar  = this.soilvar && ! is.null(dim(model[[this.vnam]]))
         if (this.soilvar){
            this.soilvar = this.soilvar && dim(model[[this.vnam]])[2] == model$nzg
         }#end if (this.soilvar)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make sure this is a soil variable, but skip soil wetness.                   #
         #---------------------------------------------------------------------------------#
         if (this.soilvar && (! this.vnam %in% "soil.wetness")){
            cat("         > ",this.desc,"...","\n")

            #----- Find the interpolated layers. ------------------------------------------#
            model[[this.vnam]] = interpol( x      = model$slz
                                         , y      = model[[this.vnam]]
                                         , xout   = obser$slz
                                         , along  = 2
                                         , is.log = FALSE
                                         , method = "monoH.FC"
                                         )#end interpol
            #------------------------------------------------------------------------------#
         }#end if (this.soilvar)
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ncompvar))
      #------------------------------------------------------------------------------------#



      #----- Copy soil-related dimensions and labels. -------------------------------------#
      model$nzg              = obser$nzg
      model$slz              = obser$slz
      model$slz.key          = obser$slz.key
      model$slz.desc         = obser$slz.desc
      model$soil.fnmean.when = obser$soil.fnmean.when
      model$soil.slz         = obser$soil.slz
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Find soil wetness with the interpolated soil moisture.                        #
      #------------------------------------------------------------------------------------#
      swater             = ifelse(is.finite(obser$soil.water),model$soil.water,NA)
      model$soil.wetness = apply(X=swater,MARGIN=2,FUN=percentil,trim=0.02)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Create daily and fortnightly averages, and also compute statistics.           #
      #------------------------------------------------------------------------------------#
      cat("       ~ Aggregating model data and computing statistics...","\n")
      dist.comp = list()
      for (v in sequence(ncompvar)){
         #------ Loop over variables. -----------------------------------------------------#
         this.compvar  = compvar[[v]]
         this.vnam     = this.compvar$vnam
         this.dmean    = paste("dmean"   ,this.vnam,sep=".")
         this.fnmean   = paste("fnmean"  ,this.vnam,sep=".")
         this.fnq025   = paste("fnq025"  ,this.vnam,sep=".")
         this.fnq975   = paste("fnq975"  ,this.vnam,sep=".")
         this.measured = paste("measured",this.vnam,sep=".")
         this.desc     = this.compvar$desc
         this.unit     = this.compvar$unit
         this.sunvar   = this.compvar$sunvar
         cat("         > ",this.desc,"...","\n")
         #---------------------------------------------------------------------------------#






         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
         #     Check whether to find u*-filter dependent variables.                        #
         #---------------------------------------------------------------------------------#
         this.ust   = paste("ust"           ,this.vnam,sep=".")
         ust.dmean  = paste("ust","dmean"   ,this.vnam,sep=".")
         ust.fnmean = paste("ust","fnmean"  ,this.vnam,sep=".")
         if (this.ust %in% names(obser)){
            cat("        > Build u*-matrix for observations...","\n")


            #----- Copy the original variable to a matrix and discard unwanted data. ------#
            model[[this.ust]] = ( matrix( data     = model[[this.vnam]]
                                        , nrow     = model$nwhen
                                        , ncol     = model$nustar
                                        , dimnames = dimnames(obser$ust.nee)
                                        )#end matrix
                                + 0. * obser[[this.ust]]
                                )#end 
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Replace gap-filled values by NA, and find the daily mean.  We replace    #
            # the original values by these, so we can eliminate model values easily.       #
            #------------------------------------------------------------------------------#
            model[[ust.dmean]]  = ( qapply( X     = model[[this.ust]]
                                          , INDEX = model$today
                                          , DIM   = 1
                                          , FUN   = mean
                                          , na.rm = TRUE
                                          )#end tapply
                                  + 0. * obser[[ust.dmean]]
                                  )#end 
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Solve fortnightly means.                                                 #
            #------------------------------------------------------------------------------#
            if (sum(is.finite(c(model[[this.ust]]))) > 0){

               #----- Save number of years. -----------------------------------------------#
               year.num = sort(unique(model$year))
               n.year   = length(year.num)
               #---------------------------------------------------------------------------#



               #----- Find mean fortnightly period. ---------------------------------------#
               fnmean.1st = qapply( X      = model[[this.ust]]
                                  , INDEX  = list(model$hr.idx,model$fortnight,model$year)
                                  , DIM    = 1
                                  , FUN    = mean
                                  , na.rm  = TRUE
                                  )#end apply
               hh = match(as.numeric(dimnames(fnmean.1st)[[1]]),hour.num            )
               ff = match(as.numeric(dimnames(fnmean.1st)[[2]]),sequence(yr.ftnight))
               yy = match(as.numeric(dimnames(fnmean.1st)[[3]]),year.num            )
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Collapse fortnightly periods.  Make sure all periods are              #
               # defined in the structure, if none of them were selected, make them        #
               #NA.                                                                        #
               #---------------------------------------------------------------------------#
               fnmean.model  = array ( data     = NA
                                     , dim      = c(nhour,yr.ftnight,n.year,obser$nustar)
                                     , dimnames = list( hour.num
                                                      , sequence(yr.ftnight)
                                                      , year.num
                                                      , model$ust.filter
                                                      )#end list
                                     )#end array
               fnmean.model[hh,ff,yy,] = fnmean.1st
               fnmean.model = apply(X=fnmean.model,MARGIN=c(2,3,4),FUN=mean)
               fnmean.model = apply(X=fnmean.model,MARGIN=c(1,3)  ,FUN=mean,na.rm=TRUE)
               fnmean.model = fnmean.model + 0. * obser[[ust.fnmean]]
               #---------------------------------------------------------------------------#

            }else{
               #----- Leave everything empty. ---------------------------------------------#
               fnmean.model = matrix( data     = NA
                                    , nrow     = yr.ftnight
                                    , ncol     = model$nustar
                                    , dimnames = list( sequence(yr.ftnight)
                                                     , model$ust.filter)
                                    )#end fnmean.obser
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#


            #----- Save the data. ---------------------------------------------------------#
            model[[ust.fnmean]] = fnmean.model
            #------------------------------------------------------------------------------#
         }#end if (this.ust %in% names(obser))
         #---------------------------------------------------------------------------------#






         #---------------------------------------------------------------------------------#
         #     Get rid of model results from times where no observation had been measured. #
         #---------------------------------------------------------------------------------#
         obs.now     = obser[[this.vnam]]
         mod.now     = model[[this.vnam]] + 0. * obser[[this.vnam]]
         res.now     = obs.now - mod.now
         dmean.model = qapply(X=mod.now,INDEX=model$today,DIM=1,FUN=mean,na.rm=TRUE)
         dmean.model = dmean.model + 0. * obser[[this.dmean]]
         #---------------------------------------------------------------------------------#




         #----- Find the number of layers. ------------------------------------------------#
         nlyrs   = ncol(obs.now)
         lyr.key = colnames(obs.now)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Solve fortnightly means and find confidence intervals for the mean.         #
         #---------------------------------------------------------------------------------#
         #----- For loop over columns (sorry, Steve...) -----------------------------------#
         fnmean.model = NULL
         fnq025.model = NULL
         fnq975.model = NULL
         for (cc in sequence(nlyrs)){
            keep  = is.finite(mod.now[,cc])

            if (any(keep)){
               #----- Find mean fortnightly period. ---------------------------------------#
               data.in    = data.frame( x         = mod.now[,cc]
                                      , fortnight = model$fortnight
                                      , year      = model$year
                                      , hour      = model$hr.idx
                                      )#end data.frame
               if (use.hier.boot){
                  fn.boot = bhier.fortnight.mean( data.in = data.in
                                                , R       = ftnight.n.boot
                                                , ci      = 0.95
                                                )#end bhier.fortnight.mean
               }else{
                  fn.boot = boot.fortnight.mean ( data.in = data.in
                                                , R       = ftnight.n.boot
                                                , ci      = 0.95
                                                )#end boot.fortnight.mean
               }#end if
               fnmean.now = fn.boot$expected
               fnq025.now = fn.boot$qlow
               fnq975.now = fn.boot$qhigh
               #---------------------------------------------------------------------------#
            }else{
               fnmean.now   = rep(x=NA,times=obser$nfnmean)
               fnq025.now   = rep(x=NA,times=obser$nfnmean)
               fnq975.now   = rep(x=NA,times=obser$nfnmean)
            }#end if
            #------------------------------------------------------------------------------#


            #------ Append bootstrap for this layer. --------------------------------------#
            fnmean.model = cbind(fnmean.model,fnmean.now)
            fnq025.model = cbind(fnq025.model,fnq025.now)
            fnq975.model = cbind(fnq975.model,fnq975.now)
            rm(fnmean.now,fnq025.now,fnq975.now)
            #------------------------------------------------------------------------------#
         }#end for (cc in sequence(nlyrs))
         #---------------------------------------------------------------------------------#




         #----- Save the data. ------------------------------------------------------------#
         model[[this.vnam  ]] = mod.now     
         model[[this.dmean ]] = dmean.model 
         model[[this.fnmean]] = fnmean.model
         model[[this.fnq025]] = fnq025.model
         model[[this.fnq975]] = fnq975.model
         #---------------------------------------------------------------------------------#






         #---------------------------------------------------------------------------------#
         #     Initialise the list to keep comparisons between model and observations.     #
         #---------------------------------------------------------------------------------#
         arr3 = array( data     = NA
                     , dim      = c(nlyrs,ndiel,nseason)
                     , dimnames = list(lyr.key,diel.key,season.key)
                     )#end matrix
         arr4 = array( data     = NA
                     , dim      = c(nlyrs,ndiel,nseason,nmoment)
                     , dimnames = list(lyr.key,diel.key,season.key,moment.key)
                     )#end array
         comp = list ( residuals   = res.now
                     , n           = arr3
                     , obs.moment  = arr4
                     , mod.moment  = arr4
                     , res.moment  = arr4
                     , bias        = arr3
                     , sigma       = arr3
                     , lsq.lnlike  = arr3
                     , mse         = arr3
                     , rmse        = arr3
                     , r.squared   = arr3
                     , fvue        = arr3
                     , sw.stat     = arr3
                     , sw.p.value  = arr3
                     , ks.stat     = arr3
                     , ks.p.value  = arr3
                     )#end list
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #   Loop over layers.                                                             #
         #---------------------------------------------------------------------------------#
         for (cc in sequence(nlyrs)){
            #------------------------------------------------------------------------------#
            #     Season block.                                                            #
            #------------------------------------------------------------------------------#
            for (ee in sequence(nseason)){
               #---------------------------------------------------------------------------#
               #    Diel block.                                                            #
               #---------------------------------------------------------------------------#
               for (dd in sequence(ndiel)){
                  #------------------------------------------------------------------------#
                  #      Select data.                                                      #
                  #------------------------------------------------------------------------#
                  if (dd %in% diel.fnmean){
                     #----- Daily means, no daytime check. --------------------------------#
                     sel     = obser$fnmean.season == ee | ee == nseason
                     obs.use = obser[[this.fnmean]][sel,cc]
                     mod.use = model[[this.fnmean]][sel,cc]
                     #---------------------------------------------------------------------#
                  }else if (dd %in% diel.dmean){
                     #----- Daily means, no daytime check. --------------------------------#
                     sel     = obser$dmean.season == ee | ee == nseason
                     obs.use = obser[[this.dmean]][sel,cc]
                     mod.use = model[[this.dmean]][sel,cc]
                     #---------------------------------------------------------------------#
                  }else{
                     #----- Skip check for sun variables during the night. ----------------#
                     e.sel   = obser$season == ee | ee == nseason
                     d.sel   = obser$diel   == dd | dd %in% diel.all.hrs
                     s.sel   = obser$highsun | ( ! this.sunvar )
                     sel     = e.sel & d.sel & s.sel
                     obs.use = obser[[this.vnam]][sel,cc]
                     mod.use = model[[this.vnam]][sel,cc]
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #      Skip statistics if no data has been selected or if everything     #
                  # is the same.                                                           #
                  #------------------------------------------------------------------------#
                  sd.obser = sd(obs.use,na.rm=TRUE)
                  #------------------------------------------------------------------------#

                  if (is.finite(sd.obser) && sd.obser > 1.0e-6){
                     #---------------------------------------------------------------------#
                     #      Find the summary of how good (or bad...) the model is.         #
                     #---------------------------------------------------------------------#
                     goodness = test.goodness(x.mod=mod.use,x.obs=obs.use)

                     comp$n         [cc,dd,ee ] = goodness$n
                     comp$obs.moment[cc,dd,ee,] = goodness$obs.moment
                     comp$mod.moment[cc,dd,ee,] = goodness$mod.moment
                     comp$res.moment[cc,dd,ee,] = goodness$res.moment
                     comp$bias      [cc,dd,ee ] = goodness$bias
                     comp$sigma     [cc,dd,ee ] = goodness$sigma
                     comp$lsq.lnlike[cc,dd,ee ] = goodness$lsq.lnlike
                     comp$mse       [cc,dd,ee ] = goodness$mse
                     comp$rmse      [cc,dd,ee ] = goodness$rmse
                     comp$r.squared [cc,dd,ee ] = goodness$r.squared
                     comp$fvue      [cc,dd,ee ] = goodness$fvue
                     comp$sw.stat   [cc,dd,ee ] = goodness$sw.statistic
                     comp$sw.p.value[cc,dd,ee ] = goodness$sw.p.value
                     comp$ks.stat   [cc,dd,ee ] = goodness$ks.statistic
                     comp$ks.p.value[cc,dd,ee ] = goodness$ks.p.value
                     #---------------------------------------------------------------------#
                  }#end if (is.finite(sd.obser) && sd.obser > 1.0e-6)
                  #------------------------------------------------------------------------#
               }#end for (dd in sequence(ndiel))
               #---------------------------------------------------------------------------#
            }#end for (ee in sequence(nseason))
            #------------------------------------------------------------------------------#
         }#end for (cc in sequence(nlyrs))
         #---------------------------------------------------------------------------------#


         #----- Save the comparison list for this variable. -------------------------------#
         dist.comp[[this.vnam]] = comp
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ncompvar))
      #------------------------------------------------------------------------------------#




      #----- Save the data and free some memory. ------------------------------------------#
      this$sim[[simul$name[s]]] = dist.comp
      this$ans[[simul$name[s]]] = model
      rm(list=c("dist.comp","model","eddy.complete","eddy.tresume"))
      #------------------------------------------------------------------------------------#
   }#end for (s in sequence(nsimul))
   #---------------------------------------------------------------------------------------#





   #----- Copy the data to the results. ---------------------------------------------------#
   eft.iata    = paste("eft",iata,sep=".")
   res.iata    = paste("res",iata,sep=".")
   assign(eft.iata,obser)
   assign(res.iata,this )
   eft[[iata]] = obser
   res[[iata]] = this
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #      Save hourly data to RData.                                                       #
   #---------------------------------------------------------------------------------------#
   rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
   cat(" + Saving hourly data to ",basename(rdata.iata),"...","\n")
   dummy = save(list=c(eft.iata,res.iata), file=rdata.iata)
   rm(eft.iata,res.iata,obser,this)
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
#      Plot the light response curve by site and season.  Each plot has one panel for each #
# simulation.                                                                              #
#------------------------------------------------------------------------------------------#
if (plot.gpp.light){
   cat(" + Find the light response curve...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over the seasons.                                                            #
   #---------------------------------------------------------------------------------------#
   for (ee in sequence(nseason)){
      cat("   - ",season.full[ee],"...","\n",sep="")

      list.obs.use  = list()
      list.obs.pred = list()
      list.mod.use  = list()
      list.mod.pred = list()
      light.xlimit  = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      light.ylimit  = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)

      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information and site observation. ---------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Select observed data sets.                                                  #
         #---------------------------------------------------------------------------------#
         if (use.dmean.light){
            #------ Keep only data with both GPP and PAR. ---------------------------------#
            f.sel   = is.finite(c(obser$par)) & is.finite(c(obser$gpp))
            gpp.now = ifelse(f.sel,c(obser$gpp),NA)
            par.now = ifelse(f.sel,c(obser$par),NA)
            #------------------------------------------------------------------------------#


            #------ Select only daytime data for this season. -----------------------------#
            e.sel   = obser$season == ee | ee == nseason
            d.sel   = obser$highsun
            sel     = e.sel & f.sel
            gpp.now = gpp.now    [sel]
            par.now = par.now    [sel]
            today   = obser$today[sel]
            #------------------------------------------------------------------------------#



            #----- Compute the daily means. -----------------------------------------------#
            dmean.gpp = tapply(X      = gpp.now
                              ,INDEX  = today
                              ,FUN    = ifenough
                              ,f      = mean
                              ,ef.min = ef.min
                              ,na.rm  = TRUE
                              )#end tapply
            dmean.par = tapply(X      = par.now
                              ,INDEX  = today
                              ,FUN    = ifenough
                              ,f      = mean
                              ,ef.min = ef.min
                              ,na.rm  = TRUE
                              )#end tapply
            #------------------------------------------------------------------------------#


            #----- Keep only the data that is useful. -------------------------------------#
            dm.sel  = is.finite(dmean.gpp) & is.finite(dmean.par)
            n.sel   = sum(dm.sel)
            obs.use = data.frame( par = dmean.par[dm.sel]
                                , gpp = dmean.gpp[dm.sel]
                                )#end data.frame
            n.fit.min = n.fit.dmean.min
            #------------------------------------------------------------------------------#
         }else{
            e.sel   = obser$season == ee | ee == nseason
            f.sel   = is.finite(c(obser$par)) & is.finite(c(obser$gpp))
            d.sel   = obser$highsun
            sel     = e.sel & f.sel & d.sel
            n.sel   = sum(sel)
            obs.use = data.frame( par = c(obser$par)[sel]
                                , gpp = c(obser$gpp)[sel]
                                )#end data.frame
            n.fit.min = n.fit.fmean.min
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fit a curve only if there are enough points.                                #
         #---------------------------------------------------------------------------------#
         if (n.sel >= n.fit.min){
            cat("     * Site :",this.longname," (N = ",n.sel,")...","\n")



            #----- Select and sort the data. ----------------------------------------------#
            if (tolower(light.method) == "nls"){
               obs.pred = nls.light.response( par.in  = obs.use$par
                                            , gpp     = obs.use$gpp
                                            , first   = c(a1=1,a2=10,a3=500)
                                            , n.boot  = light.n.boot
                                            , control = list( maxiter   = 500
                                                            , minFactor = 2^-26
                                                            )#end list
                                            )#end fit.light.response
            }else{
               obs.pred = optim.light.response( par.in  = obs.use$par
                                              , gpp     = obs.use$gpp
                                              , first   = c(a1=1,a2=10,a3=500)
                                              , n.boot  = light.n.boot
                                              , skew    = skew.optim
                                              )#end fit.light.response
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Grab the data for each simulation.                                       #
            #------------------------------------------------------------------------------#
            mod.use  = list()
            mod.pred = list()
            for (s in sequence(nsimul)){
               model = res[[iata]]$ans[[s]]

               #---------------------------------------------------------------------------#
               #     Find out when this output variable is finite and measured (or at      #
               # least shortwave radiation had been measured).                             #
               #---------------------------------------------------------------------------#
               if (use.dmean.light){

                  #------ Keep only data with both GPP and PAR. ---------------------------#
                  gpp.now = ifelse(f.sel,c(model$gpp),NA)
                  par.now = ifelse(f.sel,c(model$par),NA)
                  #------------------------------------------------------------------------#


                  #------ Select only daytime data for this season. -----------------------#
                  e.sel   = obser$season == ee | ee == nseason
                  gpp.now = gpp.now    [sel]
                  par.now = par.now    [sel]
                  today   = model$today[sel]
                  #------------------------------------------------------------------------#



                  #----- Compute the daily means. -----------------------------------------#
                  dmean.gpp = tapply(X      = gpp.now
                                    ,INDEX  = today
                                    ,FUN    = ifenough
                                    ,f      = mean
                                    ,ef.min = ef.min
                                    ,na.rm  = TRUE
                                    )#end tapply
                  dmean.par = tapply(X      = par.now
                                    ,INDEX  = today
                                    ,FUN    = ifenough
                                    ,f      = mean
                                    ,ef.min = ef.min
                                    ,na.rm  = TRUE
                                    )#end tapply
                  #------------------------------------------------------------------------#


                  mod.use[[s]] = data.frame( par = dmean.par[dm.sel]
                                           , gpp = dmean.gpp[dm.sel]
                                           )#end data.frame
               }else{
                  mod.use[[s]] = data.frame( par = c(model$par)[sel]
                                           , gpp = c(model$gpp)[sel]
                                           )#end data.frame
               }#end if
               #---------------------------------------------------------------------------#



               #----- Select and sort the data. -------------------------------------------#
               if (tolower(light.method) == "nls"){
                  mod.pred[[s]] = nls.light.response( par.in  = mod.use[[s]]$par
                                                    , gpp     = mod.use[[s]]$gpp
                                                    , first   = c(a1=1,a2=10,a3=500)
                                                    , n.boot  = light.n.boot
                                                    , control = list( maxiter   = 500
                                                                    , minFactor = 2^-26
                                                                    )#end list
                                                    )#end fit.light.response
               }else{
                  mod.pred[[s]] = optim.light.response( par.in  = mod.use[[s]]$par
                                                      , gpp     = mod.use[[s]]$gpp
                                                      , first   = c(a1=1,a2=10,a3=500)
                                                      , skew    = skew.optim
                                                      , n.boot  = light.n.boot
                                                      )#end fit.light.response
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Save range.                                                           #
               #---------------------------------------------------------------------------#
               light.xlimit[[s]] = range( c(light.xlimit[[s]],obs.use$par,mod.use[[s]]$par)
                                        , finite = TRUE
                                        )#end range
               light.ylimit[[s]] = range( c(light.ylimit[[s]],obs.use$gpp,mod.use[[s]]$gpp)
                                        , finite = TRUE
                                        )#end range
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            list.obs.use [[iata]] = obs.use
            list.obs.pred[[iata]] = obs.pred
            list.mod.use [[iata]] = mod.use
            list.mod.pred[[iata]] = mod.pred
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot simulations                                                                #
      #------------------------------------------------------------------------------------#
      xlimit = pretty.xylim(u=unlist(light.xlimit),fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=unlist(light.ylimit),fracexp=0.0,is.log=FALSE)
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         #---------------------------------------------------------------------------------#

         if (iata %in% names(list.obs.use)){
            obs.use  = list.obs.use [[iata]]
            obs.pred = list.obs.pred[[iata]]
            mod.use  = list.mod.use [[iata]]
            mod.pred = list.mod.pred[[iata]]
            n.sel    = length(obs.use$par)


            #------ Set some common features. ---------------------------------------------#
            letitre = paste("Light response curve: ",this.longname,"\n"
                           ,season.full[ee]," ( N = ",n.sel,")",sep="")
            lex     = desc.unit(desc="Incoming PAR",unit=untab$umolom2os)
            ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$kgcom2oyr)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the light response curves.                                         #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.light = out[[outform[o]]]$light$sites
               fichier = file.path(out.light,paste("gpp-light-",season.suffix[ee],"-"
                                                  ,iata,".",outform[o],sep="")
                                  )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(1,1,4,0))
               layout(mat = lo.simul$mat)
               #---------------------------------------------------------------------------#



               #----- Loop over all simulations. ------------------------------------------#
               for (s in sequence(nsimul)){
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  if (lo.simul$bottom[s]) axis.rt(side=1,las=5,off=0.1)
                  if (lo.simul$left  [s]) axis(side=2,las=1)
                  grid(col=grid.colour,lty="dotted")
                  abline(h=0,v=0,lty="solid")
                  #------ Plot the points. ------------------------------------------------#
                  points( x    = obs.use$par
                        , y    = obs.use$gpp
                        , type = "p"
                        , pch  = 16
                        , col  = obs.col
                        , cex  = 0.5
                        )#end points
                  points( x    = mod.use[[s]]$par
                        , y    = mod.use[[s]]$gpp
                        , type = "p"
                        , pch  = 16
                        , col  = ed22.col
                        , cex  = 0.5
                        )#end points
                  #------------------------------------------------------------------------#


                  #----- Plot the lines. --------------------------------------------------#
                  lines(x=obs.pred$par,y=obs.pred$gpp ,col=obs.fg,lwd=2,lty="solid" )
                  lines(x=obs.pred$par,y=obs.pred$q025,col=obs.fg,lwd=2,lty="dashed")
                  lines(x=obs.pred$par,y=obs.pred$q975,col=obs.fg,lwd=2,lty="dashed")
                  lines(x   = mod.pred[[s]]$par
                       ,y   = mod.pred[[s]]$gpp
                       ,col = ed22.fg
                       ,lwd = 2
                       ,lty = "solid"
                       )#end lines
                  lines(x   = mod.pred[[s]]$par
                       ,y   = mod.pred[[s]]$q025
                       ,col = ed22.fg
                       ,lwd = 2
                       ,lty = "dashed"
                       )#end lines
                  lines(x   = mod.pred[[s]]$par
                       ,y   = mod.pred[[s]]$q975
                       ,col = ed22.fg
                       ,lwd = 2
                       ,lty = "dashed"
                       )#end lines
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s])
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot title. ---------------------------------------------------------#
               gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
         }#end if (iata %in% names(list.obs(use)))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot default simulation.                                                        #
      #------------------------------------------------------------------------------------#
      s      = sim.default
      xlimit = pretty.xylim(u=light.xlimit[[s]],fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=light.ylimit[[s]],fracexp=0.0,is.log=FALSE)


      #------ Set some common features. ---------------------------------------------------#
      letitre = paste("Light response curve - ",season.full[ee],sep="")
      lex     = desc.unit(desc="Incoming PAR",unit=untab$umolom2os)
      ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$kgcom2oyr)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Plot the light response curves.                                               #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.light = out[[outform[o]]]$light$default
         fichier = file.path(out.light,paste("gpp-light-",season.suffix[ee]
                                            ,"-",simul$name[s],".",outform[o],sep="")
                            )#end file.path
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(1,1,4,0))
         layout(mat = lo.site$mat)
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #----- Retrieve data. ---------------------------------------------------------#
            if (iata %in% names(list.obs.use)){
               obs.use  = list.obs.use [[iata]]
               obs.pred = list.obs.pred[[iata]]
               mod.use  = list.mod.use [[iata]][[s]]
               mod.pred = list.mod.pred[[iata]][[s]]
               n.sel    = length(obs.use$par)
            }else{
               n.sel    = 0
            }#end if
            #------------------------------------------------------------------------------#



            #---- Open window. ------------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            if (lo.site$bottom[p]) axis.rt(side=1,las=5,off=0.1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline(h=0,v=0,lty="solid")
            #------------------------------------------------------------------------------#

            if (n.sel > 0){
               #------ Plot the points. ---------------------------------------------------#
               points( x    = obs.use$par
                     , y    = obs.use$gpp
                     , type = "p"
                     , pch  = 16
                     , col  = obs.col
                     , cex  = 0.5
                     )#end points
               points( x    = mod.use$par
                     , y    = mod.use$gpp
                     , type = "p"
                     , pch  = 16
                     , col  = ed22.col
                     , cex  = 0.5
                     )#end points
               #---------------------------------------------------------------------------#


               #----- Plot the lines. -----------------------------------------------------#
               lines(x=obs.pred$par,y=obs.pred$gpp ,col=obs.fg,lwd=2,lty="solid" )
               lines(x=obs.pred$par,y=obs.pred$q025,col=obs.fg,lwd=2,lty="dashed")
               lines(x=obs.pred$par,y=obs.pred$q975,col=obs.fg,lwd=2,lty="dashed")
               lines(x   = mod.pred$par
                    ,y   = mod.pred$gpp
                    ,col = ed22.fg
                    ,lwd = 2
                    ,lty = "solid"
                    )#end lines
               lines(x   = mod.pred$par
                    ,y   = mod.pred$q025
                    ,col = ed22.fg
                    ,lwd = 2
                    ,lty = "dashed"
                    )#end lines
               lines(x   = mod.pred$par
                    ,y   = mod.pred$q975
                    ,col = ed22.fg
                    ,lwd = 2
                    ,lty = "dashed"
                    )#end lines
               #---------------------------------------------------------------------------#
            }#end if (n.sel > 0)
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(sites$desc[p]," (",toupper(sites$iata[p]),")"
                            ,"\n","N = ",n.sel,sep=""),line=0.25)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.75)
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
   }#end for (ee in sequence(nseason))
   #---------------------------------------------------------------------------------------#
}#end if (plot.gpp.light)
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
#      Plot the VPD response curve by site and season.  Each plot has one panel for each   #
# simulation.                                                                              #
#------------------------------------------------------------------------------------------#
if (plot.gpp.vpdef){
   cat(" + Find the VPD response curve...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over the seasons.                                                            #
   #---------------------------------------------------------------------------------------#
   for (ee in sequence(nseason)){
      cat("   - ",season.full[ee],"...","\n",sep="")

      list.obs.use  = list()
      list.mod.use  = list()
      vpdef.xlimit  = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      vpdef.ylimit  = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)

      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information and site observation. ---------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Select observed data sets.                                                  #
         #---------------------------------------------------------------------------------#
         if (use.dmean.vpdef){
            #------ Keep only data with both GPP and PAR. ---------------------------------#
            f.sel     = is.finite(c(obser$atm.vpdef)) & is.finite(c(obser$gpp))
            gpp.now   = ifelse(f.sel,c(obser$gpp      ),NA)
            vpdef.now = ifelse(f.sel,c(obser$atm.vpdef),NA)
            #------------------------------------------------------------------------------#


            #------ Select only daytime data for this season. -----------------------------#
            e.sel     = obser$season == ee | ee == nseason
            d.sel     = obser$highsun
            sel       = e.sel & f.sel
            gpp.now   = gpp.now    [sel]
            vpdef.now = vpdef.now  [sel]
            today     = obser$today[sel]
            #------------------------------------------------------------------------------#



            #----- Compute the daily means. -----------------------------------------------#
            dmean.gpp   = tapply( X      = gpp.now
                                , INDEX  = today
                                , FUN    = ifenough
                                , f      = mean
                                , ef.min = ef.min
                                , na.rm  = TRUE
                                )#end tapply
            dmean.vpdef = tapply( X      = vpdef.now
                                , INDEX  = today
                                , FUN    = ifenough
                                , f      = mean
                                , ef.min = ef.min
                                , na.rm  = TRUE
                                )#end tapply
            #------------------------------------------------------------------------------#


            #----- Keep only the data that is useful. -------------------------------------#
            dm.sel  = is.finite(dmean.gpp) & is.finite(dmean.vpdef)
            n.sel   = sum(dm.sel)
            obs.use = data.frame( vpdef = dmean.vpdef[dm.sel]
                                , gpp   = dmean.gpp  [dm.sel]
                                )#end data.frame
            n.fit.min = n.fit.dmean.min
            #------------------------------------------------------------------------------#
         }else{
            e.sel   = obser$season == ee | ee == nseason
            f.sel   = is.finite(c(obser$atm.vpdef)) & is.finite(c(obser$gpp))
            d.sel   = obser$highsun
            sel     = e.sel & f.sel & d.sel
            n.sel   = sum(sel)
            obs.use = data.frame( vpdef = c(obser$atm.vpdef)[sel]
                                , gpp   = c(obser$gpp      )[sel]
                                )#end data.frame
            n.fit.min = n.fit.fmean.min
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fit a curve only if there are enough points.                                #
         #---------------------------------------------------------------------------------#
         if (n.sel >= n.fit.min){
            cat("     * Site :",this.longname," (N = ",n.sel,")...","\n")

            #------------------------------------------------------------------------------#
            #     Grab the data for each simulation.                                       #
            #------------------------------------------------------------------------------#
            mod.use  = list()
            mod.pred = list()
            for (s in sequence(nsimul)){
               model           = res[[iata]]$ans[[s]]
               model$atm.vpdef = ( 0.01 * vpdefil( pres = model$atm.prss * 100.
                                                 , temp = model$atm.temp + t00
                                                 , humi = model$atm.shv  * 0.001 )
                                 + 0. * obser$atm.vpdef )



               #---------------------------------------------------------------------------#
               #     Find out when this output variable is finite and measured (or at      #
               # least shortwave radiation had been measured).                             #
               #---------------------------------------------------------------------------#
               if (use.dmean.vpdef){

                  #------ Keep only data with both GPP and PAR. ---------------------------#
                  gpp.now   = ifelse(f.sel,c(model$gpp      ),NA)
                  vpdef.now = ifelse(f.sel,c(model$atm.vpdef),NA)
                  #------------------------------------------------------------------------#


                  #------ Select only daytime data for this season. -----------------------#
                  e.sel     = obser$season == ee | ee == nseason
                  gpp.now   = gpp.now    [sel]
                  vpdef.now = vpdef.now  [sel]
                  today     = model$today[sel]
                  #------------------------------------------------------------------------#



                  #----- Compute the daily means. -----------------------------------------#
                  dmean.gpp   = tapply( X      = gpp.now
                                      , INDEX  = today
                                      , FUN    = ifenough
                                      , f      = mean
                                      , ef.min = ef.min
                                      , na.rm  = TRUE
                                      )#end tapply
                  dmean.vpdef = tapply( X      = vpdef.now
                                      , INDEX  = today
                                      , FUN    = ifenough
                                      , f      = mean
                                      , ef.min = ef.min
                                      , na.rm  = TRUE
                                      )#end tapply
                  #------------------------------------------------------------------------#


                  mod.use[[s]] = data.frame( vpdef = dmean.vpdef[dm.sel]
                                           , gpp   = dmean.gpp[dm.sel]
                                           )#end data.frame
               }else{
                  mod.use[[s]] = data.frame( vpdef = c(model$atm.vpdef)[sel]
                                           , gpp   = c(model$gpp      )[sel]
                                           )#end data.frame
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Save range.                                                           #
               #---------------------------------------------------------------------------#
               vpdef.xlimit[[s]] = range( c( vpdef.xlimit[[s]]
                                           , obs.use$vpdef
                                           , mod.use[[s]]$vpdef
                                           )#end c
                                        , finite = TRUE
                                        )#end range
               vpdef.ylimit[[s]] = range( c(vpdef.ylimit[[s]],obs.use$gpp,mod.use[[s]]$gpp)
                                        , finite = TRUE
                                        )#end range
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            list.obs.use [[iata]] = obs.use
            list.mod.use [[iata]] = mod.use
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot simulations                                                                #
      #------------------------------------------------------------------------------------#
      xlimit = pretty.xylim(u=unlist(vpdef.xlimit),fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=unlist(vpdef.ylimit),fracexp=0.0,is.log=FALSE)
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         #---------------------------------------------------------------------------------#

         if (iata %in% names(list.obs.use)){
            obs.use  = list.obs.use [[iata]]
            mod.use  = list.mod.use [[iata]]
            n.sel    = length(obs.use$vpdef)


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"\n",season.full[ee]," ( N = ",n.sel,")",sep="")
            lex     = desc.unit(desc="Above-canopy VPD",unit=untab$hpa)
            ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$kgcom2oyr)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the light response curves.                                         #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.vpdef = out[[outform[o]]]$vpdef$sites
               fichier   = file.path(out.vpdef,paste("gpp-vpdef-",season.suffix[ee]
                                                    ,"-",iata,".",outform[o],sep="")
                                    )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(1,1,4,0))
               layout(mat = lo.simul$mat)
               #---------------------------------------------------------------------------#



               #----- Loop over all simulations. ------------------------------------------#
               for (s in sequence(nsimul)){
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  if (lo.simul$bottom[s]) axis(side=1)
                  if (lo.simul$left  [s]) axis(side=2,las=1)
                  grid(col=grid.colour,lty="dotted")
                  abline(h=0,v=0,lty="solid")
                  #------ Plot the points. ------------------------------------------------#
                  now.vpdef  = c(obs.use$vpdef,mod.use[[s]]$vpdef)
                  now.gpp    = c(obs.use$gpp  ,mod.use[[s]]$gpp  )
                  now.colour = c( rep(obs.col ,times=length(obs.use$vpdef))
                                , rep(ed22.col,times=length(mod.use[[s]]$vpdef))
                                )#end c
                  shuffle    = sample(length(now.vpdef))
                  points( x    = now.vpdef [shuffle]
                        , y    = now.gpp   [shuffle]
                        , type = "p"
                        , pch  = 16
                        , col  = now.colour[shuffle]
                        , cex  = 0.5
                        )#end points
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s])
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot title. ---------------------------------------------------------#
               gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
         }#end if (iata %in% names(list.obs(use)))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot default simulation.                                                        #
      #------------------------------------------------------------------------------------#
      s      = sim.default
      xlimit = pretty.xylim(u=vpdef.xlimit[[s]],fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=vpdef.ylimit[[s]],fracexp=0.0,is.log=FALSE)


      #------ Set some common features. ---------------------------------------------------#
      letitre = paste(season.full[ee],sep="")
      lex     = desc.unit(desc="Above-canopy VPD",unit=untab$hpa)
      ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$kgcom2oyr)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Plot the light response curves.                                               #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.vpdef = out[[outform[o]]]$vpdef$default
         fichier   = file.path(out.vpdef,paste("gpp-vpdef-",season.suffix[ee]
                                              ,"-",simul$name[s],".",outform[o],sep="")
                            )#end file.path
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(1,1,4,0))
         layout(mat = lo.site$mat)
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #----- Retrieve data. ---------------------------------------------------------#
            if (iata %in% names(list.obs.use)){
               obs.use  = list.obs.use [[iata]]
               mod.use  = list.mod.use [[iata]][[s]]
               n.sel    = length(obs.use$vpdef)
            }else{
               n.sel    = 0
            }#end if
            #------------------------------------------------------------------------------#



            #---- Open window. ------------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            if (lo.site$bottom[p]) axis(side=1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline(h=0,v=0,lty="solid")
            #------------------------------------------------------------------------------#




            #---- Plot GPP as a function of VPD. ------------------------------------------#
            if (n.sel > 0){
               now.vpdef  = c(obs.use$vpdef,mod.use$vpdef)
               now.gpp    = c(obs.use$gpp  ,mod.use$gpp  )
               now.colour = c( rep(obs.col ,times=length(obs.use$vpdef))
                             , rep(ed22.col,times=length(mod.use$vpdef))
                             )#end c
               shuffle    = sample(length(now.vpdef))
               points( x    = now.vpdef [shuffle]
                     , y    = now.gpp   [shuffle]
                     , type = "p"
                     , pch  = 16
                     , col  = now.colour[shuffle]
                     , cex  = 0.5
                     )#end points
               #---------------------------------------------------------------------------#
            }#end if (n.sel > 0)
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(sites$desc[p]," (",toupper(sites$iata[p]),")"
                            ,"\n","N = ",n.sel,sep=""),line=0.25)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
   }#end for (ee in sequence(nseason))
   #---------------------------------------------------------------------------------------#
}#end if (plot.gpp.vpdef)
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
#      Plot the soil moisture response curve by site and season.  Each plot has one panel  #
# for each simulation.                                                                     #
#------------------------------------------------------------------------------------------#
if (plot.gpp.wetness){
   cat(" + Find the soil moisture response curve...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over the seasons.                                                            #
   #---------------------------------------------------------------------------------------#
   for (ee in sequence(nseason)){
      cat("   - ",season.full[ee],"...","\n",sep="")

      list.obs.use   = list()
      list.mod.use   = list()
      list.cc        = list()
      wetness.xlimit = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      wetness.ylimit = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)

      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information and site observation. ---------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         ndat          = colSums(is.finite(obser$soil.wetness))
         nlyr          = length(ndat)
         slz           = obser$slz
         slz.use       = ifelse(ndat==0,-slz,slz)
         cc            = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         list.cc[[iata]] = list(cc=cc, slz = paste( "z = ",obser$slz[cc],"m", sep=""))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Select observed data sets.                                                  #
         #---------------------------------------------------------------------------------#
         if (use.dmean.wetness){
            #------ Keep only data with both GPP and PAR. ---------------------------------#
            f.sel       = is.finite(c(obser$soil.wetness[,cc])) & is.finite(c(obser$gpp))
            gpp.now     = ifelse(f.sel,c(obser$gpp      ),NA)
            wetness.now = ifelse(f.sel,c(obser$soil.wetness[,cc]),NA)
            #------------------------------------------------------------------------------#


            #------ Select only daytime data for this season. -----------------------------#
            e.sel       = obser$season == ee | ee == nseason
            d.sel       = obser$highsun
            sel         = e.sel & f.sel
            gpp.now     = gpp.now    [sel]
            wetness.now = wetness.now[sel]
            today       = obser$today[sel]
            #------------------------------------------------------------------------------#



            #----- Compute the daily means. -----------------------------------------------#
            if (length(gpp.now) > 0){
               dmean.gpp     = tapply( X      = gpp.now
                                     , INDEX  = today
                                     , FUN    = ifenough
                                     , f      = mean
                                     , ef.min = ef.min
                                     , na.rm  = TRUE
                                     )#end tapply
               dmean.wetness = tapply( X      = wetness.now
                                     , INDEX  = today
                                     , FUN    = ifenough
                                     , f      = mean
                                     , ef.min = ef.min
                                     , na.rm  = TRUE
                                     )#end tapply
               dm.sel  = is.finite(dmean.gpp) & is.finite(dmean.wetness)
               n.sel   = sum(dm.sel)
               obs.use = data.frame( wetness = dmean.wetness[dm.sel]
                                   , gpp     = dmean.gpp    [dm.sel]
                                   )#end data.frame
               #---------------------------------------------------------------------------#
            }else{
               #----- Skip this site .-----------------------------------------------------#
               n.sel = 0
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #----- Keep only the data that is useful. -------------------------------------#
            n.fit.min = n.fit.dmean.min
            #------------------------------------------------------------------------------#
         }else{
            e.sel   = obser$season == ee | ee == nseason
            f.sel   = is.finite(c(obser$soil.wetness[,cc])) & is.finite(c(obser$gpp))
            d.sel   = obser$highsun
            sel     = e.sel & f.sel & d.sel
            n.sel   = sum(sel)
            obs.use = data.frame( wetness = obser$soil.wetness[sel,cc]
                                , gpp     = c(obser$gpp      )[sel]
                                )#end data.frame
            n.fit.min = n.fit.fmean.min
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Fit a curve only if there are enough points.                                #
         #---------------------------------------------------------------------------------#
         if (n.sel >= n.fit.min){
            cat("     * Site :",this.longname," (N = ",n.sel,")...","\n")

            #------------------------------------------------------------------------------#
            #     Grab the data for each simulation.                                       #
            #------------------------------------------------------------------------------#
            mod.use  = list()
            mod.pred = list()
            for (s in sequence(nsimul)){
               model           = res[[iata]]$ans[[s]]

               #---------------------------------------------------------------------------#
               #     Find out when this output variable is finite and measured (or at      #
               # least shortwave radiation had been measured).                             #
               #---------------------------------------------------------------------------#
               if (use.dmean.wetness){

                  #------ Keep only data with both GPP and PAR. ---------------------------#
                  gpp.now     = ifelse(f.sel,c(model$gpp              ),NA)
                  wetness.now = ifelse(f.sel,c(model$soil.wetness[,cc]),NA)
                  #------------------------------------------------------------------------#


                  #------ Select only daytime data for this season. -----------------------#
                  e.sel       = obser$season == ee | ee == nseason
                  gpp.now     = gpp.now    [sel]
                  wetness.now = wetness.now[sel]
                  today       = model$today[sel]
                  #------------------------------------------------------------------------#



                  #----- Compute the daily means. -----------------------------------------#
                  dmean.gpp     = tapply( X      = gpp.now
                                        , INDEX  = today
                                        , FUN    = ifenough
                                        , f      = mean
                                        , ef.min = ef.min
                                        , na.rm  = TRUE
                                        )#end tapply
                  dmean.wetness = tapply( X      = wetness.now
                                        , INDEX  = today
                                        , FUN    = ifenough
                                        , f      = mean
                                        , ef.min = ef.min
                                        , na.rm  = TRUE
                                        )#end tapply
                  #------------------------------------------------------------------------#


                  mod.use[[s]] = data.frame( wetness = dmean.wetness[dm.sel]
                                           , gpp     = dmean.gpp    [dm.sel]
                                           )#end data.frame
               }else{
                  mod.use[[s]] = data.frame( wetness = model$soil.wetness[sel,cc]
                                           , gpp     = model$gpp         [sel,1]
                                           )#end data.frame
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Save range.                                                           #
               #---------------------------------------------------------------------------#
               wetness.xlimit[[s]] = range( c( wetness.xlimit[[s]]
                                             , obs.use$wetness
                                             , mod.use[[s]]$wetness
                                             )#end c
                                          , finite = TRUE
                                          )#end range
               wetness.ylimit[[s]] = range( c(wetness.ylimit[[s]]
                                             ,obs.use$gpp
                                             ,mod.use[[s]]$gpp
                                             )#end c
                                          , finite = TRUE
                                          )#end range
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            list.obs.use [[iata]] = obs.use
            list.mod.use [[iata]] = mod.use
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot simulations                                                                #
      #------------------------------------------------------------------------------------#
      xlimit = pretty.xylim(u=unlist(wetness.xlimit),fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=unlist(wetness.ylimit),fracexp=0.0,is.log=FALSE)
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         #---------------------------------------------------------------------------------#

         if (iata %in% names(list.obs.use)){
            obs.use  = list.obs.use [[iata]]
            mod.use  = list.mod.use [[iata]]
            cc.use   = list.cc      [[iata]]
            n.sel    = length(obs.use$wetness)


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"  -  ",cc.use$slz,"\n",season.full[ee]
                           ," ( N = ",n.sel,")",sep="")
            lex     = desc.unit(desc="Relative soil moisture"    ,unit=untab$empty    )
            ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$kgcom2oyr)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the light response curves.                                         #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.wetness = out[[outform[o]]]$wetness$sites
               fichier     = file.path(out.wetness,paste("gpp-wetness-",season.suffix[ee]
                                                        ,"-",iata,".",outform[o],sep="")
                                      )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(1,1,4,0))
               layout(mat = lo.simul$mat)
               #---------------------------------------------------------------------------#



               #----- Loop over all simulations. ------------------------------------------#
               for (s in sequence(nsimul)){
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  if (lo.simul$bottom[s]) axis(side=1)
                  if (lo.simul$left  [s]) axis(side=2,las=1)
                  grid(col=grid.colour,lty="dotted")
                  abline(h=0,v=0,lty="solid")
                  #------ Plot the points. ------------------------------------------------#
                  now.wetness  = c(obs.use$wetness,mod.use[[s]]$wetness)
                  now.gpp      = c(obs.use$gpp  ,mod.use[[s]]$gpp  )
                  now.colour   = c( rep(obs.col ,times=length(obs.use$wetness))
                                  , rep(ed22.col,times=length(mod.use[[s]]$wetness))
                                  )#end c
                  shuffle      = sample(length(now.wetness))
                  points( x    = now.wetness[shuffle]
                        , y    = now.gpp    [shuffle]
                        , type = "p"
                        , pch  = 16
                        , col  = now.colour [shuffle]
                        , cex  = 0.5
                        )#end points
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s],line=0.5)
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot title. ---------------------------------------------------------#
               gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
         }#end if (iata %in% names(list.obs(use)))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot default simulation.                                                        #
      #------------------------------------------------------------------------------------#
      s      = sim.default
      xlimit = pretty.xylim(u=wetness.xlimit[[s]],fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=wetness.ylimit[[s]],fracexp=0.0,is.log=FALSE)


      #------ Set some common features. ---------------------------------------------------#
      letitre = paste(season.full[ee],sep="")
      lex     = desc.unit(desc="Soil wetness",unit=untab$empty)
      ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$kgcom2oyr)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Plot the light response curves.                                               #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.wetness = out[[outform[o]]]$wetness$default
         fichier   = file.path(out.wetness,paste("gpp-wetness-",season.suffix[ee]
                                                ,"-",simul$name[s],".",outform[o],sep="")
                              )#end file.path
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(1,1,4,0))
         layout(mat = lo.site$mat)
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #----- Retrieve data. ---------------------------------------------------------#
            if (iata %in% names(list.obs.use)){
               obs.use  = list.obs.use [[iata]]
               mod.use  = list.mod.use [[iata]][[s]]
               n.sel    = length(obs.use$wetness)
            }else{
               n.sel    = 0
            }#end if
            cc.use      = list.cc      [[iata]]
            #------------------------------------------------------------------------------#



            #---- Open window. ------------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            if (lo.site$bottom[p]) axis(side=1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline(h=0,v=0,lty="solid")
            #------------------------------------------------------------------------------#




            #---- Plot GPP as a function of VPD. ------------------------------------------#
            if (n.sel > 0){
               now.wetness  = c(obs.use$wetness,mod.use$wetness)
               now.gpp      = c(obs.use$gpp    ,mod.use$gpp    )
               now.colour   = c( rep(obs.col ,times=length(obs.use$wetness))
                               , rep(ed22.col,times=length(mod.use$wetness))
                               )#end c
               shuffle      = sample(length(now.wetness))
               points( x    = now.wetness [shuffle]
                     , y    = now.gpp   [shuffle]
                     , type = "p"
                     , pch  = 16
                     , col  = now.colour[shuffle]
                     , cex  = 0.5
                     )#end points
               #---------------------------------------------------------------------------#
            }#end if (n.sel > 0)
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(sites$desc[p]," (",toupper(sites$iata[p]),")"
                            ,"\n",cc.use$slz,"  -  N = ",n.sel,sep=""),line=0.25)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
   }#end for (ee in sequence(nseason))
   #---------------------------------------------------------------------------------------#
}#end if (plot.gpp.wetness)
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
#      Plot the soil moisture response curve by site and season.  Each plot has one panel  #
# for each simulation.                                                                     #
#------------------------------------------------------------------------------------------#
if (plot.reco.wetness){
   cat(" + Find the Reco soil moisture response curve...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over the seasons.                                                            #
   #---------------------------------------------------------------------------------------#
   for (ee in sequence(nseason)){
      cat("   - ",season.full[ee],"...","\n",sep="")

      list.obs.use   = list()
      list.mod.use   = list()
      list.cc        = list()
      wetness.xlimit = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      wetness.ylimit = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)

      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information and site observation. ---------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         ndat          = colSums(is.finite(obser$soil.wetness))
         nlyr          = length(ndat)
         slz           = obser$slz
         slz.use       = ifelse(ndat==0,-slz,slz)
         cc            = pmin(nlyr,which.min(abs(slz.use-slz.reco.ref)))
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         list.cc[[iata]] = list(cc=cc, slz = paste( "z = ",obser$slz[cc],"m", sep=""))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Select observed data sets.                                                  #
         #---------------------------------------------------------------------------------#
         if (use.dmean.wetness){
            #------ Keep only data with both GPP and PAR. ---------------------------------#
            f.sel       = is.finite(c(obser$soil.wetness[,cc])) & is.finite(c(obser$reco))
            reco.now    = ifelse(f.sel,c(obser$reco             ),NA)
            wetness.now = ifelse(f.sel,c(obser$soil.wetness[,cc]),NA)
            #------------------------------------------------------------------------------#


            #------ Select only daytime data for this season. -----------------------------#
            e.sel       = obser$season == ee | ee == nseason
            sel         = e.sel & f.sel
            reco.now    = reco.now   [sel]
            wetness.now = wetness.now[sel]
            today       = obser$today[sel]
            #------------------------------------------------------------------------------#



            #----- Compute the daily means. -----------------------------------------------#
            if (length(reco.now) > 0){
               dmean.reco    = tapply( X      = reco.now
                                     , INDEX  = today
                                     , FUN    = ifenough
                                     , f      = mean
                                     , ef.min = ef.min
                                     , na.rm  = TRUE
                                     )#end tapply
               dmean.wetness = tapply( X      = wetness.now
                                     , INDEX  = today
                                     , FUN    = ifenough
                                     , f      = mean
                                     , ef.min = ef.min
                                     , na.rm  = TRUE
                                     )#end tapply
               dm.sel  = is.finite(dmean.reco) & is.finite(dmean.wetness)
               n.sel   = sum(dm.sel)
               obs.use = data.frame( wetness = dmean.wetness[dm.sel]
                                   , reco    = dmean.reco   [dm.sel]
                                   )#end data.frame
               #---------------------------------------------------------------------------#
            }else{
               #----- Skip this site .-----------------------------------------------------#
               n.sel = 0
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #----- Keep only the data that is useful. -------------------------------------#
            n.fit.min = n.fit.dmean.min
            #------------------------------------------------------------------------------#
         }else{
            e.sel   = obser$season == ee | ee == nseason
            f.sel   = is.finite(c(obser$soil.wetness[,cc])) & is.finite(c(obser$reco))
            sel     = e.sel & f.sel
            n.sel   = sum(sel)
            obs.use = data.frame( wetness = obser$soil.wetness[sel,cc]
                                , reco    = c(obser$reco     )[sel]
                                )#end data.frame
            n.fit.min = n.fit.fmean.min
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Fit a curve only if there are enough points.                                #
         #---------------------------------------------------------------------------------#
         if (n.sel >= n.fit.min){
            cat("     * Site :",this.longname," (N = ",n.sel,")...","\n")

            #------------------------------------------------------------------------------#
            #     Grab the data for each simulation.                                       #
            #------------------------------------------------------------------------------#
            mod.use  = list()
            mod.pred = list()
            for (s in sequence(nsimul)){
               model           = res[[iata]]$ans[[s]]

               #---------------------------------------------------------------------------#
               #     Find out when this output variable is finite and measured (or at      #
               # least shortwave radiation had been measured).                             #
               #---------------------------------------------------------------------------#
               if (use.dmean.wetness){

                  #------ Keep only data with both GPP and PAR. ---------------------------#
                  reco.now    = ifelse(f.sel,c(model$reco             ),NA)
                  wetness.now = ifelse(f.sel,c(model$soil.wetness[,cc]),NA)
                  #------------------------------------------------------------------------#


                  #------ Select only daytime data for this season. -----------------------#
                  e.sel       = obser$season == ee | ee == nseason
                  reco.now    = reco.now   [sel]
                  wetness.now = wetness.now[sel]
                  today       = model$today[sel]
                  #------------------------------------------------------------------------#



                  #----- Compute the daily means. -----------------------------------------#
                  dmean.reco    = tapply( X      = reco.now
                                        , INDEX  = today
                                        , FUN    = ifenough
                                        , f      = mean
                                        , ef.min = ef.min
                                        , na.rm  = TRUE
                                        )#end tapply
                  dmean.wetness = tapply( X      = wetness.now
                                        , INDEX  = today
                                        , FUN    = ifenough
                                        , f      = mean
                                        , ef.min = ef.min
                                        , na.rm  = TRUE
                                        )#end tapply
                  #------------------------------------------------------------------------#


                  mod.use[[s]] = data.frame( wetness = dmean.wetness[dm.sel]
                                           , reco    = dmean.reco   [dm.sel]
                                           )#end data.frame
               }else{
                  mod.use[[s]] = data.frame( wetness = model$soil.wetness[sel,cc]
                                           , reco    = model$reco        [sel,1]
                                           )#end data.frame
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Save range.                                                           #
               #---------------------------------------------------------------------------#
               wetness.xlimit[[s]] = range( c( wetness.xlimit[[s]]
                                             , obs.use$wetness
                                             , mod.use[[s]]$wetness
                                             )#end c
                                          , finite = TRUE
                                          )#end range
               wetness.ylimit[[s]] = range( c(wetness.ylimit[[s]]
                                             ,obs.use$reco
                                             ,mod.use[[s]]$reco
                                             )#end c
                                          , finite = TRUE
                                          )#end range
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            list.obs.use [[iata]] = obs.use
            list.mod.use [[iata]] = mod.use
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot simulations                                                                #
      #------------------------------------------------------------------------------------#
      xlimit = pretty.xylim(u=unlist(wetness.xlimit),fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=unlist(wetness.ylimit),fracexp=0.0,is.log=FALSE)
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         #---------------------------------------------------------------------------------#

         if (iata %in% names(list.obs.use)){
            obs.use  = list.obs.use [[iata]]
            mod.use  = list.mod.use [[iata]]
            cc.use   = list.cc      [[iata]]
            n.sel    = length(obs.use$wetness)


            #------ Set some common features. ---------------------------------------------#
            letitre = paste(this.longname,"  -  ",cc.use$slz,"\n",season.full[ee]
                           ," ( N = ",n.sel,")",sep="")
            lex     = desc.unit(desc="Relative soil moisture",unit=untab$empty    )
            ley     = desc.unit(desc="Ecosystem respiration" ,unit=untab$kgcom2oyr)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the light response curves.                                         #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.wetness = out[[outform[o]]]$wetness$sites
               fichier     = file.path(out.wetness,paste("reco-wetness-",season.suffix[ee]
                                                        ,"-",iata,".",outform[o],sep="")
                                      )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(1,1,4,0))
               layout(mat = lo.simul$mat)
               #---------------------------------------------------------------------------#



               #----- Loop over all simulations. ------------------------------------------#
               for (s in sequence(nsimul)){
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  if (lo.simul$bottom[s]) axis(side=1)
                  if (lo.simul$left  [s]) axis(side=2,las=1)
                  grid(col=grid.colour,lty="dotted")
                  abline(h=0,v=0,lty="solid")
                  #------ Plot the points. ------------------------------------------------#
                  now.wetness  = c(obs.use$wetness,mod.use[[s]]$wetness)
                  now.reco     = c(obs.use$reco   ,mod.use[[s]]$reco   )
                  now.colour   = c( rep(obs.col ,times=length(obs.use$wetness))
                                  , rep(ed22.col,times=length(mod.use[[s]]$wetness))
                                  )#end c
                  shuffle      = sample(length(now.wetness))
                  points( x    = now.wetness[shuffle]
                        , y    = now.reco   [shuffle]
                        , type = "p"
                        , pch  = 16
                        , col  = now.colour [shuffle]
                        , cex  = 0.5
                        )#end points
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s],line=0.5)
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot title. ---------------------------------------------------------#
               gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
         }#end if (iata %in% names(list.obs(use)))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Plot default simulation.                                                        #
      #------------------------------------------------------------------------------------#
      s      = sim.default
      xlimit = pretty.xylim(u=wetness.xlimit[[s]],fracexp=0.0,is.log=FALSE)
      ylimit = pretty.xylim(u=wetness.ylimit[[s]],fracexp=0.0,is.log=FALSE)


      #------ Set some common features. ---------------------------------------------------#
      letitre = paste(season.full[ee],sep="")
      lex     = desc.unit(desc="Soil wetness"         ,unit=untab$empty    )
      ley     = desc.unit(desc="Ecosystem respiration",unit=untab$kgcom2oyr)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Plot the light response curves.                                               #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.wetness = out[[outform[o]]]$wetness$default
         fichier   = file.path(out.wetness,paste("reco-wetness-",season.suffix[ee]
                                                ,"-",simul$name[s],".",outform[o],sep="")
                              )#end file.path
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(1,1,4,0))
         layout(mat = lo.site$mat)
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #----- Retrieve data. ---------------------------------------------------------#
            if (iata %in% names(list.obs.use)){
               obs.use  = list.obs.use [[iata]]
               mod.use  = list.mod.use [[iata]][[s]]
               n.sel    = length(obs.use$wetness)
            }else{
               n.sel    = 0
            }#end if
            cc.use      = list.cc      [[iata]]
            #------------------------------------------------------------------------------#



            #---- Open window. ------------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            if (lo.site$bottom[p]) axis(side=1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline(h=0,v=0,lty="solid")
            #------------------------------------------------------------------------------#




            #---- Plot GPP as a function of VPD. ------------------------------------------#
            if (n.sel > 0){
               now.wetness  = c(obs.use$wetness,mod.use$wetness)
               now.reco     = c(obs.use$reco   ,mod.use$reco   )
               now.colour   = c( rep(obs.col ,times=length(obs.use$wetness))
                               , rep(ed22.col,times=length(mod.use$wetness))
                               )#end c
               shuffle      = sample(length(now.wetness))
               points( x    = now.wetness[shuffle]
                     , y    = now.reco   [shuffle]
                     , type = "p"
                     , pch  = 16
                     , col  = now.colour[shuffle]
                     , cex  = 0.5
                     )#end points
               #---------------------------------------------------------------------------#
            }#end if (n.sel > 0)
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(sites$desc[p]," (",toupper(sites$iata[p]),")"
                            ,"\n",cc.use$slz,"  -  N = ",n.sel,sep=""),line=0.25)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=3.5)
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
   }#end for (ee in sequence(nseason))
   #---------------------------------------------------------------------------------------#
}#end if (plot.wetness)
#------------------------------------------------------------------------------------------#
























#------------------------------------------------------------------------------------------#
#      Plot the time series of fortnightly means.                                          #
#------------------------------------------------------------------------------------------#
if (plot.ust.ftnight){
   cat(" + Mapping u*-filter dependent variables onto the u*-filter map...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.dmean      = paste(      "dmean" ,this.vnam,sep=".")
      this.fnmean     = paste(      "fnmean",this.vnam,sep=".")
      this.fnq025     = paste(      "fnq025",this.vnam,sep=".")
      this.fnq975     = paste(      "fnq975",this.vnam,sep=".")
      this.ust        = paste("ust"         ,this.vnam,sep=".")
      this.ust.fnmean = paste("ust","fnmean",this.vnam,sep=".")
      this.ust.fnq025 = paste("ust","fnq025",this.vnam,sep=".")
      this.ust.fnq975 = paste("ust","fnq975",this.vnam,sep=".")
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      is.ust          = this.compvar$ustvar
      this.soilvar    = this.compvar$soilvar
      this.sunvar     = this.compvar$sunvar
      this.col        = this.compvar$col
      this.fg         = this.compvar$fg
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high

      #------------------------------------------------------------------------------------#
      #      Check whether this is a u*-dependent variable.                                #
      #------------------------------------------------------------------------------------#
      if (is.ust){

         cat("   - ",this.desc,"...","\n")
         ustfn.obser    = list()
         ustfn.model    = list()
         sites.ylimit   = list()
         simul.ylimit   = mapply(FUN=numeric,length=0*sequence(nsimul),SIMPLIFY=FALSE)
         sites.zlimit   = list()
         sites.blimit   = list()
         simul.blimit   = mapply(FUN=numeric,length=0*sequence(nsimul),SIMPLIFY=FALSE)

         #---------------------------------------------------------------------------------#
         #     First loop: grab data.                                                      #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information and observed data. ---------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            obser         = eft[[iata]]
            cat("     * Site :",this.longname,"...","\n")
            #------------------------------------------------------------------------------#






            #------------------------------------------------------------------------------#
            #     Grab fortnightly means for this site, and find out whether there is any- #
            # thing to plot.                                                               #
            #------------------------------------------------------------------------------#
            #----- If this is NEE, fix the u* scale, I forgot to change it to umol/m2/s. --#
            if (this.vnam %in% "nee"){ mult = kgCyr.2.umols }else{ mult = 1.0 }
            obs.now       = obser[[this.vnam]] * mult
            fnmean.obser  = split( x = obser[[this.ust.fnmean]] * mult
                                 , f = row(obser[[this.ust.fnmean]])
                                 )#end split
            #------------------------------------------------------------------------------#
            cnt           = obser$highsun | (! this.sunvar)
            ndat          = sum(is.finite(obs.now[cnt,]))
            #------------------------------------------------------------------------------#



            #----- Find out what to append. -----------------------------------------------#
            ust.title = paste(this.longname," (",toupper(iata),") - N = ",ndat,sep="")
            #------------------------------------------------------------------------------#



            #------ Plot only if the data are finite. -------------------------------------#
            if (sum(ndat) > 0){
               #----- Use only the subset of data. ----------------------------------------#
               col.use               = (sequence(obser$nustar) %% (ust.row.skip+1)) == 1
               col.use[obser$nustar] = TRUE


               #----- Save the data. ------------------------------------------------------#
               ustfn.obser  [[iata]] = list( when        = obser$ust.fnmean.when
                                           , ust         = obser$ust.ustmin
                                           , ust.default = obser$ust.filter[obser$ubest]
                                           , ust.altern  = obser$ust.filter[obser$ualt ]
                                           , expected    = obser[[this.ust.fnmean]]
                                           , ust.title   = ust.title
                                           , col.use     = col.use
                                           )#end if
               sites.zlimit [[iata]] = range(unlist(fnmean.obser),finite=TRUE)
               sites.ylimit [[iata]] = range(obser$ust.filter)
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Loop over simulations, and grab the expected value and confidence    #
               # interval.                                                                 #
               #---------------------------------------------------------------------------#
               ust.best  = list()
               bias.best = list()
               for (s in sequence(nsimul)){

                  #----- Grab mean fortnightly period. ------------------------------------#
                  model        = res[[iata]]$ans[[s]]
                  fnmean.model = split( x = model[[this.ust.fnmean]]
                                      , f = row(model[[this.ust.fnmean]])
                                      )#end split
                  #------------------------------------------------------------------------#



                  #------ Find the u* that makes the best match. --------------------------#
                  ust.best [[s]] = mapply( FUN      = closest
                                         , x        = fnmean.model
                                         , A        = fnmean.obser
                                         , MoreArgs = list(y=obser$ust.filter)
                                         )#end mapply
                  bias.best[[s]] = mapply( FUN      = closest.bias
                                         , x        = fnmean.model
                                         , A        = fnmean.obser
                                         )#end mapply
                  #------------------------------------------------------------------------#



                  #----- Update range. ----------------------------------------------------#
                  simul.ylimit[[   s]] = c(simul.ylimit[[s]]   ,obser$ust.filter)
                  sites.blimit[[iata]] = c(sites.blimit[[iata]],bias.best[[s]]  )
                  simul.blimit[[   s]] = c(simul.blimit[[s]]   ,bias.best[[s]]  )
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               ustfn.model[[iata]] = list( when      = model$fnmean.when
                                         , ust.best  = ust.best
                                         , bias.best = bias.best
                                         )#end list
               #---------------------------------------------------------------------------#
            }#end if (sum(ndat) > 0)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Second loop: the plots.                                                     #
         #---------------------------------------------------------------------------------#
         cat("     * Plotting by scenario...","\n")
         loop.iata = which(sites$iata %in% names(ustfn.obser))
         for (p in loop.iata){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            obs.when        = ustfn.obser[[iata]]$when
            obs.ust         = ustfn.obser[[iata]]$ust
            obs.ust.default = ustfn.obser[[iata]]$ust.default
            obs.ust.altern  = ustfn.obser[[iata]]$ust.altern
            obs.expected    = ustfn.obser[[iata]]$expected
            ust.title       = ustfn.obser[[iata]]$ust.title
            col.use         = ustfn.obser[[iata]]$col.use
            mod.when        = ustfn.model[[iata]]$when
            ust.best        = ustfn.model[[iata]]$ust.best
            bias.best       = ustfn.model[[iata]]$bias.best
            ylimit          = pretty.xylim(u=sites.ylimit[[iata]],fracexp=0.0,is.log=FALSE)
            if (this.vnam %in% c("nee","nep")){
               zspan = c(-1.04,1.04) * quantile(abs(sites.zlimit[[iata]])
                                               ,prob=0.95,na.rm=TRUE)
            }else{
               zspan = quantile(sites.zlimit[[iata]],prob=c(0.025,0.975),na.rm=TRUE)
               zspan = zspan+c(-0.04,0.04)*diff(zspan)
            }#end if
            zlimit   = pretty.xylim(u=zspan,fracexp=0.0,is.log=FALSE)
            #------------------------------------------------------------------------------#



            #------ Set some common features. ---------------------------------------------#
            letitre = ust.title
            ley     = desc.unit(desc="u* filter",unit=untab$mos)
            lex     = "Month"
            lacle   = desc.unit(desc=NULL,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------ Find which levels to plot the u-* labels. -----------------------------#
            ust.at = pretty(x=ylimit)
            if (this.vnam %in% c("nee","nep")){
               z.pal     = two.palettes(x=zlimit,n=ncol.ust,white=1
                                       ,low=hue.low,high=hue.high,zero="grey85")
               z.at      = z.pal$breaks
               z.colours = z.pal$colours
            }else{
               z.at      = pretty(x=zlimit,n=ncol.ust)
               z.colours = cscheme(n=length(z.at)-1)
            }#end if
            #------------------------------------------------------------------------------#



            #------ Find axes limits and size of the dots. --------------------------------#
            if (is.null(bias.max.std)){
               bias.max   = max(abs(sites.blimit[[iata]]))
               bias.max   = ifelse(is.na(bias.max),1,bias.max)
            }else if (this.vnam %in% "nee"){
               bias.max   = bias.max.std * kgCyr.2.umols
            }else{
               bias.max   = bias.max.std
            }#end if
            bias.scale = pretty(x=c(0,bias.max))
            bias.max   = rev(bias.scale)[1]
            bias.min   = bias.scale[2]
            cex.slope  = (ust.cex.max-ust.cex.min) / (bias.max - bias.min)
            leg.bias   = unique(c(rev(-bias.scale),bias.scale))
            leg.cex    = ust.cex.min + pmax(0,(abs(leg.bias)-bias.min)*cex.slope)
            leg.pch    = ifelse(abs(leg.bias) < bias.min,23,ifelse( leg.bias > 0,24,25))
            n.leg      = length(leg.bias)
            leg.title  = desc.unit(desc="Model bias",unit=this.unit)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Replicate the observations for each simulation, and set up the          #
            # additional settings.                                                         #
            #------------------------------------------------------------------------------#
            x.bg           = list()
            y.bg           = list()
            z.bg           = list()
            dy.bg          = list()
            x.axis.options = list()
            y.axis.options = list()
            sub.options    = list()
            plot.after     = list()
            for (s in sequence(nsimul)){
               bias.now = bias.best[[s]]

               #----- Find the symbol and the size of the simbol. -------------------------#
               cex.now  = ust.cex.min + pmax( 0,(abs(bias.now)-bias.min)*cex.slope
                                            , na.rm = TRUE
                                            )#end pmax
               pch.now  = ifelse( is.na(bias.now) | abs(bias.now) < bias.min
                                , 23
                                , ifelse( bias.now > 0, 24, 25)
                                )#end ifelse
               #---------------------------------------------------------------------------#


               x.bg          [[s]] = obs.when    [,col.use]
               y.bg          [[s]] = obs.ust     [,col.use]
               z.bg          [[s]] = obs.expected[,col.use]
               z.bg          [[s]] = pmax(zlimit[1]+tiny.num
                                         ,pmin(zlimit[2]-tiny.num,z.bg[[s]]))
               dy.bg         [[s]] = median(diff(y.bg[[s]][1,]))
               x.axis.options[[s]] = list( side     = 1
                                         , at       = fnmean.at
                                         , labels   = fnmean.labels
                                         , cex.axis = 0.8*cex.ptsz
                                         )#end list
               y.axis.options[[s]] = list(side=2,at=ust.at,las=1)
               sub.options   [[s]] = list(main=simul$desc[s],line=1.0)
               plot.after    [[s]] = list( abline = list( v    = fnmean.at
                                                        , h    = ust.at
                                                        , col  = grid.colour
                                                        , lty  = "dotted"
                                                        )#end list
                                         , abline = list( h    = obs.ust.default
                                                        , col  = col.ust.default
                                                        , lwd  = 2.0
                                                        , lty  = "dotdash"
                                                        )#end list
                                         , abline = list( h    = obs.ust.altern
                                                        , col  = col.ust.altern
                                                        , lwd  = 2.0
                                                        , lty  = "longdash"
                                                        )#end list
                                         , points = list( x    = mod.when
                                                        , y    = ust.best[[s]]
                                                        , type = "o"
                                                        , pch  = pch.now
                                                        , cex  = cex.now
                                                        , col  = foreground
                                                        , bg   = grey.mg
                                                        , lwd  = 1.25
                                                        , xpd  = TRUE
                                                        )#end list
                                         )#end list
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the fortnightly means.                                             #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ust.ftnight$sites[[this.vnam]]
               fichier = file.path(out.now,paste("ustbest-",this.vnam,"-",iata,"."
                                  ,outform[o],sep="")
                                  )#end file.path
               if (outform[o] == "x11"){
                  X11(width=tsize$width,height=tsize$height,pointsize=col.use)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=tsize$width*depth,height=tsize$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=tsize$width,height=tsize$height
                            ,pointsize=ptsz,paper=tsize$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=tsize$width,height=tsize$height
                     ,pointsize=ptsz,paper=tsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user) 
               image.map( x              = x.bg
                        , y              = y.bg
                        , z              = z.bg
                        , dy             = dy.bg
                        , xlim           = range(fnmean.at)
                        , col            = z.colours
                        , levels         = z.at
                        , na.col         = "transparent"
                        , x.axis.options = x.axis.options
                        , y.axis.options = y.axis.options
                        , sub.options    = sub.options
                        , main.title     = list( main     = letitre
                                               , xlab     = lex
                                               , ylab     = ley
                                               , cex.main = cex.main
                                               )#end list
                        , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                        , legend.options = list( x       = "bottom"
                                               , inset  = 0
                                               , legend = leg.bias
                                               , pch    = leg.pch
                                               , pt.cex = leg.cex
                                               , ncol   = n.leg
                                               , title  = leg.title
                                               , bty    = "n"
                                               , xpd    = TRUE
                                               )#end list
                        , key.log        = FALSE
                        , key.vertical   = TRUE
                        , matrix.plot    = TRUE
                        , plot.after     = plot.after
                        , f.key          = ust.key.frac
                        , f.leg          = ust.leg.frac
                        , off.xlab       = twothirds * ust.leg.frac
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
         #     Second loop: the plots.                                                     #
         #---------------------------------------------------------------------------------#
         cat("     * Plotting default scenario...","\n")
         s      = sim.default
         ylimit = pretty.xylim(u=simul.ylimit[[s]]   ,fracexp=0.0,is.log=FALSE)
         if (this.vnam %in% c("nee","nep")){
            zspan = c(-1.04,1.04) * quantile(abs(unlist(sites.zlimit))
                                            ,prob=0.95,na.rm=TRUE)
         }else{
            zspan = quantile(unlist(sites.zlimit),c(0.025,0.975),finite=TRUE)
            zspan = zspan+c(-0.04,0.04)*diff(zspan)
         }#end if
         zlimit   = pretty.xylim(u=zspan,fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#


         #------ Find axes limits and size of the dots. -----------------------------------#
         if (is.null(bias.max.std)){
            bias.max   = max(abs(simul.blimit[[s]]))
            bias.max   = ifelse(is.na(bias.max),1,bias.max)
         }else if (this.vnam %in% "nee"){
            bias.max   = bias.max.std * kgCyr.2.umols
         }else{
            bias.max   = bias.max.std
         }#end if
         bias.max   = ifelse(is.na(bias.max),1,bias.max)
         bias.scale = pretty(x=c(0,bias.max))
         bias.max   = rev(bias.scale)[1]
         bias.min   = bias.scale[2]
         cex.slope  = (ust.cex.max-ust.cex.min) / (bias.max - bias.min)
         leg.bias   = unique(c(rev(-bias.scale),bias.scale))
         leg.cex    = ust.cex.min + pmax(0,(abs(leg.bias)-bias.min)*cex.slope)
         leg.pch    = ifelse(abs(leg.bias) < bias.min, 23, ifelse( leg.bias > 0, 24, 25))
         n.leg      = length(leg.bias)
         leg.title  = desc.unit(desc="Model bias",unit=this.unit)
         #---------------------------------------------------------------------------------#



         #------ Set some common features. ------------------------------------------------#
         letitre = this.desc
         ley     = desc.unit(desc="u* filter",unit=untab$mos)
         lex     = "Month"
         lacle   = desc.unit(desc=NULL,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #------ Find which levels to plot the u-* labels. --------------------------------#
         ust.at = pretty(x=ylimit)
         if (this.vnam %in% c("nee","nep")){
            z.pal     = two.palettes(x=zlimit,n=ncol.ust,white=1
                                    ,low=hue.low,high=hue.high,zero="grey85")
            z.at      = z.pal$breaks
            z.colours = z.pal$colours
         }else{
            z.at      = pretty(x=zlimit,n=ncol.ust)
            z.colours = cscheme(n=length(z.at)-1)
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Plot the fortnightly means.                                                #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ust.ftnight$default
            fichier = file.path(out.now,paste("ustbest-",this.vnam,"-",simul$name[s],"."
                                             ,outform[o],sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=tsize$width,height=tsize$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=tsize$width*depth,height=tsize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=tsize$width,height=tsize$height
                         ,pointsize=ptsz,paper=tsize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=tsize$width,height=tsize$height
                  ,pointsize=ptsz,paper=tsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Initialise list with all dat.                                           #
            #------------------------------------------------------------------------------#
            x.bg           = list()
            y.bg           = list()
            z.bg           = list()
            dy.bg          = list()
            x.axis.options = list()
            y.axis.options = list()
            sub.options    = list()
            plot.after     = list()
            #------------------------------------------------------------------------------#



            #----- Loop over all sites. ---------------------------------------------------#
            for (p in sequence(nsites)){



               #----- Get the basic information. ------------------------------------------#
               iata          = sites$iata[p]
               this.longname = sites$desc[p]
               plotit        = iata %in% names(ustfn.obser)
               if (plotit){
                  obs.when        = ustfn.obser[[iata]]$when
                  obs.ust         = ustfn.obser[[iata]]$ust
                  obs.ust.default = ustfn.obser[[iata]]$ust.default
                  obs.ust.altern  = ustfn.obser[[iata]]$ust.altern
                  obs.expected    = ustfn.obser[[iata]]$expected
                  ust.title       = ustfn.obser[[iata]]$ust.title
                  col.use         = ustfn.obser[[iata]]$col.use
                  mod.when        = ustfn.model[[iata]]$when
                  ust.best        = ustfn.model[[iata]]$ust.best [[s]]
                  bias.best       = ustfn.model[[iata]]$bias.best[[s]]

                  #----- Find the symbol and the size of the simbol. ----------------------#
                  cex.now  = ust.cex.min + pmax( 0,(abs(bias.best)-bias.min)*cex.slope
                                               , na.rm = TRUE
                                               )#end pmax
                  pch.now  = ifelse( is.na(bias.best) | abs(bias.best) < bias.min
                                   , 23
                                   , ifelse( bias.best > 0, 24, 25)
                                   )#end ifelse
                  #------------------------------------------------------------------------#

                  x.bg          [[p]] = obs.when    [,col.use]
                  y.bg          [[p]] = obs.ust     [,col.use]
                  dy.bg         [[p]] = median(diff(y.bg[[p]][1,]))
                  z.bg          [[p]] = obs.expected[,col.use]
                  z.bg          [[p]] = pmax(zlimit[1]+tiny.num
                                         ,pmin(zlimit[2]-tiny.num,z.bg[[p]]))

                  x.axis.options[[p]] = list( side     = 1
                                            , at       = fnmean.at
                                            , labels   = fnmean.labels
                                            , cex.axis = 0.8*cex.ptsz
                                         )#end list
                  y.axis.options[[p]] = list(side=2,at=ust.at,las=1)
                  sub.options   [[p]] = list(main=ust.title,line=1.0)
                  plot.after    [[p]] = list( abline = list( v    = fnmean.at
                                                           , h    = ust.at
                                                           , col  = grid.colour
                                                           , lty  = "dotted"
                                                           )#end list
                                            , abline = list( h    = obs.ust.default
                                                           , col  = col.ust.default
                                                           , lwd  = 2.0
                                                           , lty  = "dotdash"
                                                           )#end list
                                            , abline = list( h    = obs.ust.altern
                                                           , col  = col.ust.altern
                                                           , lwd  = 2.0
                                                           , lty  = "longdash"
                                                           )#end list
                                            , points = list( x    = mod.when
                                                           , y    = ust.best
                                                           , pch  = pch.now
                                                           , cex  = cex.now
                                                           , col  = foreground
                                                           , bg   = grey.mg
                                                           , type = "o"
                                                           , lwd  = 1.25
                                                           , xpd  = TRUE
                                                           )#end list
                                            )#end list
               }else{
                  ust.title = paste(this.longname," (",toupper(iata),") - N = 0",sep="")
                  col.use = (sequence(eft[[iata]]$nustar) %% (ust.row.skip+1)) == 1
                  col.use[eft[[iata]]$nustar] = TRUE
                  x.bg[[p]] = eft[[iata]]$ust.fnmean.when[,col.use]
                  y.bg[[p]] = eft[[iata]]$ust.ustmin     [,col.use]
                  z.bg[[p]] = eft[[iata]]$ust.ustmin     [,col.use] * NA
                  x.axis.options[[p]] = list(side=1,at=fnmean.at,labels=fnmean.labels)
                  y.axis.options[[p]] = list(side=2,at=ust.at,las=1)
                  sub.options   [[p]] = list(main=ust.title,line=1.0)
                  plot.after    [[p]] = list( abline = list( v = fnmean.at
                                                           , h = ust.at
                                                           , col = grid.colour
                                                           , lty = "dotted"
                                                           )#end list
                                            )#end list
               }#end if (plotit)
               #---------------------------------------------------------------------------#
            }#end for p in sequence(nplaces)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Split the device into multiple panels.                                   #
            #------------------------------------------------------------------------------#
            par(par.user)
            image.map( x              = x.bg
                     , y              = y.bg
                     , z              = z.bg
                     , dy             = dy.bg
                     , xlim           = range(fnmean.at)
                     , col            = z.colours
                     , levels         = z.at
                     , na.col         = "transparent"
                     , x.axis.options = x.axis.options
                     , y.axis.options = y.axis.options
                     , sub.options    = sub.options
                     , main.title     = list( main     = letitre
                                            , xlab     = lex
                                            , ylab     = ley
                                            , cex.main = cex.main
                                            )#end list
                     , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                     , legend.options = list( x       = "bottom"
                                            , inset  = 0
                                            , legend = leg.bias
                                            , pch    = leg.pch
                                            , pt.cex = leg.cex
                                            , ncol   = n.leg
                                            , title  = leg.title
                                            , bty    = "n"
                                            , xpd    = TRUE
                                            )#end list
                     , key.log        = FALSE
                     , key.vertical   = TRUE
                     , matrix.plot    = TRUE
                     , plot.after     = plot.after
                     , f.key          = ust.key.frac
                     , f.leg          = ust.leg.frac
                     , off.xlab       = twothirds * ust.leg.frac
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
      }#end if (is.ust)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ust.ftnight)
#------------------------------------------------------------------------------------------#
























#------------------------------------------------------------------------------------------#
#      Plot the time series of fortnightly means.                                          #
#------------------------------------------------------------------------------------------#
if (plot.ust.bias){
   cat(" + Mapping u*-filter differences onto the u*-filter map...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.dmean      = paste(      "dmean" ,this.vnam,sep=".")
      this.fnmean     = paste(      "fnmean",this.vnam,sep=".")
      this.fnq025     = paste(      "fnq025",this.vnam,sep=".")
      this.fnq975     = paste(      "fnq975",this.vnam,sep=".")
      this.ust        = paste("ust"         ,this.vnam,sep=".")
      this.ust.fnmean = paste("ust","fnmean",this.vnam,sep=".")
      this.ust.fnq025 = paste("ust","fnq025",this.vnam,sep=".")
      this.ust.fnq975 = paste("ust","fnq975",this.vnam,sep=".")
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      is.ust          = this.compvar$ustvar
      this.soilvar    = this.compvar$soilvar
      this.sunvar     = this.compvar$sunvar
      this.col        = this.compvar$col
      this.fg         = this.compvar$fg
      cscheme         = get(this.compvar$cscheme.mean)
      hue.low         = this.compvar$hue.low
      hue.high        = this.compvar$hue.high

      #------------------------------------------------------------------------------------#
      #      Check whether this is a u*-dependent variable.                                #
      #------------------------------------------------------------------------------------#
      if (is.ust){

         cat("   - ",this.desc,"...","\n")
         ustfn.obser    = list()
         ustfn.model    = list()
         sites.ylimit   = list()
         simul.ylimit   = mapply(FUN=numeric,length=0*sequence(nsimul),SIMPLIFY=FALSE)
         sites.zlimit   = list()
         simul.zlimit   = mapply(FUN=numeric,length=0*sequence(nsimul),SIMPLIFY=FALSE)

         #---------------------------------------------------------------------------------#
         #     First loop: grab data.                                                      #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information and observed data. ---------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            obser         = eft[[iata]]
            cat("     * Site :",this.longname,"...","\n")
            #------------------------------------------------------------------------------#






            #------------------------------------------------------------------------------#
            #     Grab fortnightly means for this site, and find out whether there is any- #
            # thing to plot.                                                               #
            #------------------------------------------------------------------------------#
            #----- If this is NEE, fix the u* scale, I forgot to change it to umol/m2/s. --#
            if (this.vnam %in% "nee"){ mult = kgCyr.2.umols }else{ mult = 1.0 }
            obs.now       = obser[[this.vnam]] * mult
            fnmean.obser  = split( x = obser[[this.ust.fnmean]] * mult
                                 , f = row(obser[[this.ust.fnmean]])
                                 )#end split
            #------------------------------------------------------------------------------#
            cnt           = obser$highsun | (! this.sunvar)
            ndat          = sum(is.finite(obs.now[cnt,]))
            #------------------------------------------------------------------------------#



            #----- Find out what to append. -----------------------------------------------#
            ust.title = paste(this.longname," (",toupper(iata),") - N = ",ndat,sep="")
            #------------------------------------------------------------------------------#



            #------ Plot only if the data are finite. -------------------------------------#
            if (sum(ndat) > 0){
               #----- Use only the subset of data. ----------------------------------------#
               col.use               = (sequence(obser$nustar) %% (ust.row.skip+1)) == 1
               col.use[obser$nustar] = TRUE


               #----- Save the data. ------------------------------------------------------#
               ustfn.obser  [[iata]] = list( when        = obser$ust.fnmean.when
                                           , ust         = obser$ust.ustmin
                                           , ust.default = obser$ust.filter[obser$ubest]
                                           , ust.altern  = obser$ust.filter[obser$ualt ]
                                           , expected    = obser[[this.ust.fnmean]]
                                           , ust.title   = ust.title
                                           , col.use     = col.use
                                           )#end if
               sites.ylimit [[iata]] = range(obser$ust.filter)
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Loop over simulations, and grab the expected value and confidence    #
               # interval.                                                                 #
               #---------------------------------------------------------------------------#
               ust.bias  = list()
               for (s in sequence(nsimul)){

                  #----- Grab mean fortnightly period and find the bias. ------------------#
                  model         = res[[iata]]$ans[[s]]
                  ust.bias[[s]] = model[[this.ust.fnmean]] - obser[[this.ust.fnmean]]
                  #------------------------------------------------------------------------#



                  #----- Update range. ----------------------------------------------------#
                  simul.ylimit[[   s]] = c(simul.ylimit[[s]]   ,obser$ust.filter)
                  sites.zlimit[[iata]] = c(sites.zlimit[[iata]],ust.bias[[s]]  )
                  simul.zlimit[[   s]] = c(simul.zlimit[[s]]   ,ust.bias[[s]]  )
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               ustfn.model[[iata]] = list( when      = model$fnmean.when
                                         , ust.bias  = ust.bias
                                         )#end list
               #---------------------------------------------------------------------------#
            }#end if (sum(ndat) > 0)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Second loop: the plots.                                                     #
         #---------------------------------------------------------------------------------#
         cat("     * Plotting by scenario...","\n")
         loop.iata = which(sites$iata %in% names(ustfn.obser))
         for (p in loop.iata){
            #----- Get the basic information. ---------------------------------------------#
            iata            = sites$iata[p]
            this.longname   = sites$desc[p]
            obs.when        = ustfn.obser[[iata]]$when
            obs.ust         = ustfn.obser[[iata]]$ust
            obs.ust.default = ustfn.obser[[iata]]$ust.default
            obs.ust.altern  = ustfn.obser[[iata]]$ust.altern
            obs.expected    = ustfn.obser[[iata]]$expected
            ust.title       = ustfn.obser[[iata]]$ust.title
            col.use         = ustfn.obser[[iata]]$col.use
            mod.when        = ustfn.model[[iata]]$when
            ust.bias        = ustfn.model[[iata]]$ust.bias
            ylimit          = pretty.xylim(u=sites.ylimit[[iata]],fracexp=0.0,is.log=FALSE)
            #------------------------------------------------------------------------------#


            #------ Find axes limits and size of the dots. --------------------------------#
            if (is.null(bias.max.std)){
               bias.max   = max(abs(sites.zlimit[[iata]]))
               bias.max   = ifelse(is.na(bias.max),1,bias.max)
            }else if (this.vnam %in% "nee"){
               bias.max   = bias.max.std * kgCyr.2.umols
            }else{
               bias.max   = bias.max.std
            }#end if
            zlimit = pretty.xylim(u=c(-bias.max,bias.max),fracexp=0.0,is.log=FALSE)
            #------------------------------------------------------------------------------#



            #------ Set some common features. ---------------------------------------------#
            letitre = ust.title
            ley     = desc.unit(desc="u* filter",unit=untab$mos)
            lex     = "Month"
            lacle   = desc.unit(desc=NULL,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------ Find which levels to plot the u-* labels. -----------------------------#
            ust.at    = pretty(x=ylimit)
            z.pal     = two.palettes(x=zlimit,n=ncol.ust,white=1
                                    ,low=hue.low,high=hue.high,zero="grey85")
            z.at      = z.pal$breaks
            z.colours = z.pal$colours
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Replicate the observations for each simulation, and set up the          #
            # additional settings.                                                         #
            #------------------------------------------------------------------------------#
            x.bg           = list()
            y.bg           = list()
            z.bg           = list()
            dy.bg          = list()
            x.axis.options = list()
            y.axis.options = list()
            sub.options    = list()
            plot.after     = list()
            for (s in sequence(nsimul)){


               x.bg          [[s]] = obs.when     [,col.use]
               y.bg          [[s]] = obs.ust      [,col.use]
               z.bg          [[s]] = ust.bias[[s]][,col.use]
               z.bg          [[s]] = pmax(zlimit[1]+tiny.num
                                         ,pmin(zlimit[2]-tiny.num,z.bg[[s]]))
               dy.bg         [[s]] = median(diff(y.bg[[s]][1,]))
               x.axis.options[[s]] = list( side     = 1
                                         , at       = fnmean.at
                                         , labels   = fnmean.labels
                                         , cex.axis = 0.8*cex.ptsz
                                         )#end list
               y.axis.options[[s]] = list(side=2,at=ust.at,las=1)
               sub.options   [[s]] = list(main=simul$desc[s],line=1.0)
               plot.after    [[s]] = list( abline = list( v    = fnmean.at
                                                        , h    = ust.at
                                                        , col  = grid.colour
                                                        , lty  = "dotted"
                                                        )#end list
                                         , abline = list( h    = obs.ust.default
                                                        , col  = col.ust.default
                                                        , lwd  = 2.0
                                                        , lty  = "dotdash"
                                                        )#end list
                                         , abline = list( h    = obs.ust.altern
                                                        , col  = col.ust.altern
                                                        , lwd  = 2.0
                                                        , lty  = "longdash"
                                                        )#end list
                                         )#end list
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the fortnightly means.                                             #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ust.bias$sites[[this.vnam]]
               fichier = file.path(out.now,paste("ustbias-",this.vnam,"-",iata,"."
                                  ,outform[o],sep="")
                                  )#end file.path
               if (outform[o] == "x11"){
                  X11(width=tsize$width,height=tsize$height,pointsize=col.use)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=tsize$width*depth,height=tsize$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=tsize$width,height=tsize$height
                            ,pointsize=ptsz,paper=tsize$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=tsize$width,height=tsize$height
                     ,pointsize=ptsz,paper=tsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user) 
               image.map( x              = x.bg
                        , y              = y.bg
                        , z              = z.bg
                        , dy             = dy.bg
                        , xlim           = range(fnmean.at)
                        , col            = z.colours
                        , levels         = z.at
                        , na.col         = "transparent"
                        , x.axis.options = x.axis.options
                        , y.axis.options = y.axis.options
                        , sub.options    = sub.options
                        , main.title     = list( main      = letitre
                                               , xlab      = lex
                                               , ylab      = ley
                                               , cex.main  = cex.main
                                               , line.xlab = 4.0
                                               )#end list
                        , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                        , key.log        = FALSE
                        , key.vertical   = TRUE
                        , matrix.plot    = TRUE
                        , plot.after     = plot.after
                        , f.key          = ust.key.frac
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
         #     Second loop: the plots.                                                     #
         #---------------------------------------------------------------------------------#
         cat("     * Plotting default scenario...","\n")
         s      = sim.default
         ylimit = pretty.xylim(u=simul.ylimit[[s]]   ,fracexp=0.0,is.log=FALSE)
         #------ Find axes limits and size of the dots. -----------------------------------#
         if (is.null(bias.max.std)){
            bias.max   = max(abs(unlist(sites.zlimit)))
            bias.max   = ifelse(is.na(bias.max),1,bias.max)
         }else if (this.vnam %in% "nee"){
            bias.max   = bias.max.std * kgCyr.2.umols
         }else{
            bias.max   = bias.max.std
         }#end if
         zlimit = pretty.xylim(u=c(-bias.max,bias.max),fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#



         #------ Set some common features. ------------------------------------------------#
         letitre = this.desc
         ley     = desc.unit(desc="u* filter",unit=untab$mos)
         lex     = "Month"
         lacle   = desc.unit(desc=NULL,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #------ Find which levels to plot the u-* labels. --------------------------------#
         ust.at    = pretty(x=ylimit)
         z.pal     = two.palettes(x=zlimit,n=ncol.ust,white=1
                                 ,low=hue.low,high=hue.high,zero="grey85")
         z.at      = z.pal$breaks
         z.colours = z.pal$colours
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Plot the fortnightly means.                                                #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$ust.bias$default
            fichier = file.path(out.now,paste("ustbias-",this.vnam,"-",simul$name[s],"."
                                             ,outform[o],sep="")
                               )#end file.path
            if (outform[o] == "x11"){
               X11(width=tsize$width,height=tsize$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=tsize$width*depth,height=tsize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=tsize$width,height=tsize$height
                         ,pointsize=ptsz,paper=tsize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=tsize$width,height=tsize$height
                  ,pointsize=ptsz,paper=tsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Initialise list with all dat.                                           #
            #------------------------------------------------------------------------------#
            x.bg           = list()
            y.bg           = list()
            z.bg           = list()
            dy.bg          = list()
            x.axis.options = list()
            y.axis.options = list()
            sub.options    = list()
            plot.after     = list()
            #------------------------------------------------------------------------------#



            #----- Loop over all sites. ---------------------------------------------------#
            for (p in sequence(nsites)){



               #----- Get the basic information. ------------------------------------------#
               iata          = sites$iata[p]
               this.longname = sites$desc[p]
               plotit        = iata %in% names(ustfn.obser)
               if (plotit){
                  obs.when        = ustfn.obser[[iata]]$when
                  obs.ust         = ustfn.obser[[iata]]$ust
                  obs.ust.default = ustfn.obser[[iata]]$ust.default
                  obs.ust.altern  = ustfn.obser[[iata]]$ust.altern
                  obs.expected    = ustfn.obser[[iata]]$expected
                  ust.title       = ustfn.obser[[iata]]$ust.title
                  col.use         = ustfn.obser[[iata]]$col.use
                  mod.when        = ustfn.model[[iata]]$when
                  ust.bias        = ustfn.model[[iata]]$ust.bias[[s]]

                  x.bg          [[p]] = obs.when    [,col.use]
                  y.bg          [[p]] = obs.ust     [,col.use]
                  dy.bg         [[p]] = median(diff(y.bg[[p]][1,]))
                  z.bg          [[p]] = ust.bias    [,col.use]
                  z.bg          [[p]] = pmax(zlimit[1]+tiny.num
                                         ,pmin(zlimit[2]-tiny.num,z.bg[[p]]))

                  x.axis.options[[p]] = list( side     = 1
                                            , at       = fnmean.at
                                            , labels   = fnmean.labels
                                            , cex.axis = 0.8*cex.ptsz
                                         )#end list
                  y.axis.options[[p]] = list(side=2,at=ust.at,las=1)
                  sub.options   [[p]] = list(main=ust.title,line=1.0)
                  plot.after    [[p]] = list( abline = list( v    = fnmean.at
                                                           , h    = ust.at
                                                           , col  = grid.colour
                                                           , lty  = "dotted"
                                                           )#end list
                                            , abline = list( h    = obs.ust.default
                                                           , col  = col.ust.default
                                                           , lwd  = 2.0
                                                           , lty  = "dotdash"
                                                           )#end list
                                            , abline = list( h    = obs.ust.altern
                                                           , col  = col.ust.altern
                                                           , lwd  = 2.0
                                                           , lty  = "longdash"
                                                           )#end list
                                            )#end list
               }else{
                  ust.title = paste(this.longname," (",toupper(iata),") - N = 0",sep="")
                  col.use = (sequence(eft[[iata]]$nustar) %% (ust.row.skip+1)) == 1
                  col.use[eft[[iata]]$nustar] = TRUE
                  x.bg[[p]] = eft[[iata]]$ust.fnmean.when[,col.use]
                  y.bg[[p]] = eft[[iata]]$ust.ustmin     [,col.use]
                  z.bg[[p]] = eft[[iata]]$ust.ustmin     [,col.use] * NA
                  x.axis.options[[p]] = list(side=1,at=fnmean.at,labels=fnmean.labels)
                  y.axis.options[[p]] = list(side=2,at=ust.at,las=1)
                  sub.options   [[p]] = list(main=ust.title,line=1.0)
                  plot.after    [[p]] = list( abline = list( v = fnmean.at
                                                           , h = ust.at
                                                           , col = grid.colour
                                                           , lty = "dotted"
                                                           )#end list
                                            )#end list
               }#end if (plotit)
               #---------------------------------------------------------------------------#
            }#end for p in sequence(nplaces)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Split the device into multiple panels.                                   #
            #------------------------------------------------------------------------------#
            par(par.user)
            image.map( x              = x.bg
                     , y              = y.bg
                     , z              = z.bg
                     , dy             = dy.bg
                     , xlim           = range(fnmean.at)
                     , col            = z.colours
                     , levels         = z.at
                     , na.col         = "transparent"
                     , x.axis.options = x.axis.options
                     , y.axis.options = y.axis.options
                     , sub.options    = sub.options
                     , main.title     = list( main      = letitre
                                            , xlab      = lex
                                            , ylab      = ley
                                            , cex.main  = cex.main
                                            , line.xlab = 4.0
                                            )#end list
                     , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                     , key.log        = FALSE
                     , key.vertical   = TRUE
                     , matrix.plot    = TRUE
                     , plot.after     = plot.after
                     , f.key          = ust.key.frac
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
      }#end if (is.ust)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ust.diff)
#------------------------------------------------------------------------------------------#



































#------------------------------------------------------------------------------------------#
#      Plot the time series of fortnightly means.                                          #
#------------------------------------------------------------------------------------------#
if (plot.ts.ftnight){
   cat(" + Find the fortnightly means...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.dmean      = paste("dmean"       ,this.vnam,sep=".")
      this.fnmean     = paste("fnmean"      ,this.vnam,sep=".")
      this.fnq025     = paste("fnq025"      ,this.vnam,sep=".")
      this.fnq975     = paste("fnq975"      ,this.vnam,sep=".")
      this.measured   = paste("measured"    ,this.vnam,sep=".")
      this.alt        = paste("alt"         ,this.vnam,sep=".")
      this.alt.fnmean = paste("alt","fnmean",this.vnam,sep=".")
      this.alt.fnq025 = paste("alt","fnq025",this.vnam,sep=".")
      this.alt.fnq975 = paste("alt","fnq975",this.vnam,sep=".")
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      this.soilvar    = this.compvar$soilvar
      this.sunvar     = this.compvar$sunvar
      this.col        = this.compvar$col
      this.fg         = this.compvar$fg


      cat("   - ",this.desc,"...","\n")
      tsfn.obser     = list()
      tsfn.alter     = list()
      tsfn.model     = list()
      sites.ylimit   = list()
      simul.ylimit   = mapply(FUN=numeric,length=0*sequence(nsimul),SIMPLIFY=FALSE)

      #------------------------------------------------------------------------------------#
      #     First loop: grab data.                                                         #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information and observed data. ------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         cat("     * Site :",this.longname,"...","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Grab fortnightly means for this site, and find out whether there is any-    #
         # thing to plot.                                                                  #
         #---------------------------------------------------------------------------------#
         obs.now       = obser[[this.vnam  ]]
         fnmean.obser  = obser[[this.fnmean]]
         fnq025.obser  = obser[[this.fnq025]]
         fnq975.obser  = obser[[this.fnq975]]
         cnt           = obser$highsun | (! this.sunvar)
         ndat          = colSums(is.finite(obs.now[cnt,,drop=FALSE]))
         nlyr          = ncol(fnmean.obser)
         #---------------------------------------------------------------------------------#



         #----- Find out what to append. --------------------------------------------------#
         if (this.soilvar){
            lyr.key   = obser$slz.key
            lyr.desc  = paste(this.desc,obser$slz.desc,sep="")
            lyr.title = paste(this.longname," (",toupper(iata),")",obser$slz.desc,sep="")
         }else{
            lyr.key   = ""
            lyr.desc  = this.desc
            lyr.title = paste(this.longname," (",toupper(iata),")",sep="")
         }#end if
         #---------------------------------------------------------------------------------#



         #------ Plot only if the data are finite. ----------------------------------------#
         if (sum(ndat) > 0){
            tsfn.obser  [[iata]] = list( when      = obser$fnmean.when
                                       , expected  = fnmean.obser
                                       , q025      = fnq025.obser
                                       , q975      = fnq975.obser
                                       , ndat      = ndat
                                       , nlyr      = nlyr
                                       , lyr.key   = lyr.key
                                       , lyr.desc  = lyr.desc
                                       , lyr.title = lyr.title
                                       , slz       = obser$slz
                                       )#end if
            sites.ylimit[[iata]] = range(c(fnmean.obser,fnq025.obser,fnq975.obser)
                                        ,finite=TRUE)



            #------------------------------------------------------------------------------#
            #------------------------------------------------------------------------------#
            if (this.alt %in% names(obser)){
               alt.now              = obser[[this.alt       ]]
               fnmean.alter         = obser[[this.alt.fnmean]]
               fnq025.alter         = obser[[this.alt.fnq025]]
               fnq975.alter         = obser[[this.alt.fnq975]]
               tsfn.alter  [[iata]] = list( when      = obser$fnmean.when
                                          , expected  = fnmean.alter
                                          , q025      = fnq025.alter
                                          , q975      = fnq975.alter
                                          )#end if
               sites.ylimit[[iata]] = range(c(sites.ylimit[[iata]]
                                             ,fnmean.alter,fnq025.alter,fnq975.alter)
                                           ,finite=TRUE)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            mod.expected = list()
            mod.q025     = list()
            mod.q975     = list()
            for (s in sequence(nsimul)){
               #----- Grab mean fortnightly period. ---------------------------------------#
               model             = res[[iata]]$ans[[s]]
               mod.expected[[s]] = model[[this.fnmean]]
               mod.q025    [[s]] = model[[this.fnq025]]
               mod.q975    [[s]] = model[[this.fnq975]]
               #---------------------------------------------------------------------------#


               #----- Update range. -------------------------------------------------------#
               simul.ylimit[[s]]    = range( c( simul.ylimit[[s]]
                                              , fnmean.obser
                                              , fnq025.obser
                                              , fnq975.obser
                                              , mod.expected[[s]]
                                              , mod.q025    [[s]]
                                              , mod.q975    [[s]]
                                              )#end c
                                           , finite = TRUE
                                           )#end range
               sites.ylimit[[iata]] = range( c( sites.ylimit[[iata]]
                                              , mod.expected[[   s]]
                                              , mod.q025    [[   s]]
                                              , mod.q975    [[   s]]
                                              )#end c
                                           , finite = TRUE
                                           )#end range
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            tsfn.model[[iata]] = list(expected=mod.expected,q025=mod.q025,q975=mod.q975)
            #------------------------------------------------------------------------------#
         }#end if (sum(ndat) > 0)
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Second loop: the plots.                                                        #
      #------------------------------------------------------------------------------------#
      cat("     * Plotting by scenario...","\n")
      loop.iata = which(sites$iata %in% names(tsfn.obser))
      for (p in loop.iata){
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         drya          = as.numeric(chron(paste(sites$drya[p],fnmean.year,sep="/")))
         dryz          = as.numeric(chron(paste(sites$dryz[p],fnmean.year,sep="/")))
         nlyr          = tsfn.obser[[iata]]$nlyr
         ylimit        = pretty.xylim(u=sites.ylimit[[iata]],fracexp=0.0,is.log=FALSE)
         #---------------------------------------------------------------------------------#





         #------ Make two periods in case Jan 1st is dry season. --------------------------#
         if (drya > dryz){
            drya = c(as.numeric(chron(paste(01,01,fnmean.year+00,sep="/")),drya))
            dryz = c(dryz,as.numeric(chron(paste(01,01,fnmean.year+01,sep="/"))))
         }#end if
         #---------------------------------------------------------------------------------#



         #------ Grab data for this layer. ------------------------------------------------#
         obs.when      = tsfn.obser[[iata]]$when
         obs.expected  = tsfn.obser[[iata]]$expected
         obs.q025      = tsfn.obser[[iata]]$q025
         obs.q975      = tsfn.obser[[iata]]$q975
         if (iata %in% names(tsfn.alter)){
            alt.when      = tsfn.alter[[iata]]$when
            alt.expected  = tsfn.alter[[iata]]$expected
            alt.q025      = tsfn.alter[[iata]]$q025
            alt.q975      = tsfn.alter[[iata]]$q975
         }else{
            alt.when      = NA + tsfn.obser[[iata]]$when
            alt.expected  = NA + tsfn.obser[[iata]]$expected
            alt.q025      = NA + tsfn.obser[[iata]]$q025
            alt.q975      = NA + tsfn.obser[[iata]]$q975
         }#end if
         mod.expected  = tsfn.model[[iata]]$expected
         mod.q025      = tsfn.model[[iata]]$q025
         mod.q975      = tsfn.model[[iata]]$q975
         ndat          = tsfn.obser[[iata]]$ndat
         lyr.key       = tsfn.obser[[iata]]$lyr.key
         lyr.desc      = tsfn.obser[[iata]]$lyr.desc
         lyr.title     = tsfn.obser[[iata]]$lyr.title
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over layers.                                                          #
         #---------------------------------------------------------------------------------#
         for (cc in sequence(nlyr)){
            #------ Set some common features. ---------------------------------------------#
            letitre = paste(lyr.title[cc]," N = ",ndat[cc],sep="")
            ley     = desc.unit(desc=lyr.desc[cc],unit=this.unit)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot the fortnightly means.                                             #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ts.ftnight$sites[[this.vnam]]
               fichier = file.path(out.now,paste("fnmean-",this.vnam,"-",iata,lyr.key[cc]
                                                ,".",outform[o],sep="")
                                  )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(0,1,4,0))
               layout(mat = lo.simul$mat)
               #---------------------------------------------------------------------------#



               #----- Loop over all simulations. ------------------------------------------#
               for (s in sequence(nsimul)){
                  #------ Grab model data for this simulation and column. -----------------#
                  #------------------------------------------------------------------------#



                  #------ Open sub-plot. --------------------------------------------------#
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=fnmean.limit,ylim=ylimit)
                  if (lo.simul$bottom[s]){
                     axis(side=1,at=fnmean.at,labels=fnmean.labels,cex.axis=cex.ptsz*0.75)
                  }#end if
                  if (lo.simul$left  [s]){
                     axis(side=2,las=1)
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Plot the dry season.                                              #
                  #------------------------------------------------------------------------#
                  if (show.dryseason){
                     xlimit  = as.numeric(fnmean.limit)
                     xleft   = pretty.xylim(u=xlimit,fracexp=-0.04)[1]
                     xright  = pretty.xylim(u=xlimit,fracexp=+0.04)[2]
                     ybottom = pretty.xylim(u=ylimit,fracexp=-0.04)[1]
                     ytop    = pretty.xylim(u=ylimit,fracexp=+0.04)[2]
                     xleft   = ifelse(drya <= xlimit[1],xleft ,drya)
                     xright  = ifelse(dryz >= xlimit[2],xright,dryz)
                     rect( xleft   = xleft
                         , ybottom = ybottom
                         , xright  = xright
                         , ytop    = ytop
                         , col     = col.dryseason
                         , border  = NA
                         , density = -1
                         )#end rect
                  }#end if
                  abline(v=fnmean.at,h=axTicks(2),col=grid.colour,lty="dotted")
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Plot the confidence bands (make sure that NA periods are          #
                  # properly skipped).                                                     #
                  #------------------------------------------------------------------------#
                  #----- Find the limits. -------------------------------------------------#
                  x025      = c(rbind(fnmean.when-dfnmean.when,fnmean.when+dfnmean.when))
                  alt.x025  = c(rbind(fnmean.when-dfnmean.when,fnmean.when+dfnmean.when))
                  obs.y025  = c(rbind(obs.q025     [,cc],obs.q025     [,cc]))
                  obs.y975  = c(rbind(obs.q975     [,cc],obs.q975     [,cc]))
                  alt.y025  = c(rbind(alt.q025     [,cc],alt.q025     [,cc]))
                  alt.y975  = c(rbind(alt.q975     [,cc],alt.q975     [,cc]))
                  mod.y025  = c(rbind(mod.q025[[s]][,cc],mod.q025[[s]][,cc]))
                  mod.y975  = c(rbind(mod.q975[[s]][,cc],mod.q975[[s]][,cc]))
                  keep      = ! (is.na(obs.y025) | is.na(obs.y975))
                  iblck     = cumsum(!keep)
                  kalt      = ! (is.na(alt.y025) | is.na(alt.y975))
                  ibalt     = cumsum(!kalt)
                  #----- Split polygons into finite blocks. -------------------------------#
                  x025      = split (x=x025     [keep],f=iblck[keep])
                  alt.x025  = split (x=alt.x025 [kalt],f=ibalt[kalt])
                  obs.y025  = split (x=obs.y025 [keep],f=iblck[keep])
                  mod.y025  = split (x=mod.y025 [keep],f=iblck[keep])
                  alt.y025  = split (x=alt.y025 [kalt],f=ibalt[kalt])
                  x975      = lapply(X=x025                                   ,FUN=rev)
                  alt.x975  = lapply(X=alt.x025                               ,FUN=rev)
                  obs.y975  = lapply(X=split (x=obs.y975 [keep],f=iblck[keep]),FUN=rev)
                  mod.y975  = lapply(X=split (x=mod.y975 [keep],f=iblck[keep]),FUN=rev)
                  alt.y975  = lapply(X=split (x=alt.y975 [kalt],f=ibalt[kalt]),FUN=rev)
                  #----- Plot polygons. ---------------------------------------------------#
                  npoly     = length(x025)
                  npalt     = length(alt.x025)
                  for (y in sequence(npoly)){
                     epolygon( x       = c(    x025[[y]],    x975[[y]])
                             , y       = c(obs.y025[[y]],obs.y975[[y]])
                             , col     = obs.col
                             , angle   = 45
                             , density = 20
                             )#end epolygon
                  }#end for
                  for (y in sequence(npalt)){
                     epolygon( x       = c(alt.x025[[y]],alt.x975[[y]])
                             , y       = c(alt.y025[[y]],alt.y975[[y]])
                             , col     = alt.col
                             , angle   = -45
                             , density =  20
                             )#end epolygon
                  }#end for
                  for (y in sequence(npoly)){
                     epolygon( x       = c(    x025[[y]],    x975[[y]])
                             , y       = c(mod.y025[[y]],mod.y975[[y]])
                             , col     = ed22.col
                             , angle   = +90
                             , density =  20
                             )#end epolygon
                  }#end poly
                  #------------------------------------------------------------------------#



                  #------ Plot the expected values. ---------------------------------------#
                  points( x    = fnmean.when
                        , y    = obs.expected[,cc]
                        , type = "o"
                        , pch  = 16
                        , lwd  = 2
                        , col  = obs.fg
                        )#end points
                  points( x    = fnmean.when
                        , y    = alt.expected[,cc]
                        , type = "o"
                        , pch  = 16
                        , lwd  = 2
                        , col  = alt.fg
                        )#end points
                  points( x    = fnmean.when
                        , y    = mod.expected[[s]][,cc]
                        , type = "o"
                        , pch  = 16
                        , lwd  = 2
                        , col  = ed22.fg
                        )#end points
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s],line=1,cex.main=0.8*cex.ptsz)
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot title. ---------------------------------------------------------#
               gtitle(main=letitre,ylab=ley,line.ylab=2.5)
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
         }#end for (cc in sequence(nlyr))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#









      #------------------------------------------------------------------------------------#
      #     Second loop: the plots.                                                        #
      #------------------------------------------------------------------------------------#
      cat("     * Plotting default scenario...","\n")
      s      = sim.default
      ylimit = pretty.xylim(u=simul.ylimit[[s]],fracexp=0.0,is.log=FALSE)



      #------ Set some common features. ---------------------------------------------------#
      letitre = "Fortnightly means"
      ley     = desc.unit(desc=this.desc,unit=this.unit)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Plot the fortnightly means.                                                   #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$ts.ftnight$default
         fichier = file.path(out.now,paste("fnmean-",this.vnam,"-",simul$name[s],"."
                                          ,outform[o],sep="")
                            )#end file.path
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(0,1,4,0))
         layout(mat = lo.site$mat)
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            drya          = as.numeric(chron(paste(sites$drya,fnmean.year,sep="/")))
            dryz          = as.numeric(chron(paste(sites$dryz,fnmean.year,sep="/")))
            plotit        = iata %in% names(tsfn.obser)
            if (plotit){
               ndat           = tsfn.obser[[iata]]$ndat
               nlyr           = tsfn.obser[[iata]]$nlyr
               slz            = tsfn.obser[[iata]]$slz
               slz.use        = ifelse(ndat==0,-slz,slz)
               cc             = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
               lyr.title      = tsfn.obser[[iata]]$lyr.title    [ cc]
               ndat           = tsfn.obser[[iata]]$ndat         [ cc]
               obs.expected   = tsfn.obser[[iata]]$expected     [,cc]
               obs.q025       = tsfn.obser[[iata]]$q025         [,cc]
               obs.q975       = tsfn.obser[[iata]]$q975         [,cc]
               if (iata %in% names(tsfn.alter)){
                  alt.expected   = tsfn.alter[[iata]]$expected     [,cc]
                  alt.q025       = tsfn.alter[[iata]]$q025         [,cc]
                  alt.q975       = tsfn.alter[[iata]]$q975         [,cc]
               }else{
                  alt.expected   = NA + tsfn.obser[[iata]]$expected     [,cc]
                  alt.q025       = NA + tsfn.obser[[iata]]$q025         [,cc]
                  alt.q975       = NA + tsfn.obser[[iata]]$q975         [,cc]
               }#end if
               mod.expected   = tsfn.model[[iata]]$expected[[s]][,cc]
               mod.q025       = tsfn.model[[iata]]$q025    [[s]][,cc]
               mod.q975       = tsfn.model[[iata]]$q975    [[s]][,cc]
            }else{
               ndat  = 0
               obser = eft[[iata]]
               cc    = 1

               #----- Find out what to append. --------------------------------------------#
               if (this.soilvar){
                  lyr.key   = obser$slz.key[cc]
                  lyr.desc  = paste(this.desc,obser$slz.desc[cc],sep="")
                  lyr.title = paste(this.longname,"(",toupper(iata),")"
                                   ,obser$slz.desc[cc],sep="")
               }else{
                  lyr.key   = ""
                  lyr.desc  = this.desc
                  lyr.title = paste(this.longname,"(",toupper(iata),")",sep="")
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------ Make two periods in case Jan 1st is dry season. -----------------------#
            drya          = as.numeric(chron(paste(sites$drya[p],fnmean.year,sep="/")))
            dryz          = as.numeric(chron(paste(sites$dryz[p],fnmean.year,sep="/")))
            if (drya > dryz){
               drya = c(as.numeric(chron(paste(01,01,fnmean.year+00,sep="/")),drya))
               dryz = c(dryz,as.numeric(chron(paste(01,01,fnmean.year+01,sep="/"))))
            }#end if
            #------------------------------------------------------------------------------#



            #------ Open sub-plot. --------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=fnmean.limit,ylim=ylimit)
            if (lo.site$bottom[p]){
               axis(side=1,at=fnmean.at,labels=fnmean.labels,cex.axis=cex.ptsz*0.75)
            }#end if
            if (lo.site$left  [p]){
               axis(side=2,las=1)
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot the dry season.                                                    #
            #------------------------------------------------------------------------------#
            if (show.dryseason){
               xlimit  = as.numeric(fnmean.limit)
               xleft   = pretty.xylim(u=xlimit,fracexp=-0.04)[1]
               xright  = pretty.xylim(u=xlimit,fracexp=+0.04)[2]
               ybottom = pretty.xylim(u=ylimit,fracexp=-0.04)[1]
               ytop    = pretty.xylim(u=ylimit,fracexp=+0.04)[2]
               xleft   = ifelse(drya <= xlimit[1],xleft ,drya)
               xright  = ifelse(dryz >= xlimit[2],xright,dryz)
               rect( xleft      = xleft
                   , ybottom    = ybottom
                   , xright     = xright
                   , ytop       = ytop
                   , col        = col.dryseason
                   , border     = NA
                   , density    = -1
                   )#end rect
            }#end if
            abline(v=fnmean.at,h=axTicks(2),col=grid.colour,lty="dotted")
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the confidence bands (make sure that NA periods are                #
            # properly skipped).                                                           #
            #------------------------------------------------------------------------------#
            if (plotit){

               #----- Find the limits. ----------------------------------------------------#
               x025      = c(rbind(fnmean.when-dfnmean.when,fnmean.when+dfnmean.when))
               alt.x025  = c(rbind(fnmean.when-dfnmean.when,fnmean.when+dfnmean.when))
               obs.y025  = c(rbind(obs.q025,obs.q025))
               obs.y975  = c(rbind(obs.q975,obs.q975))
               alt.y025  = c(rbind(alt.q025,alt.q025))
               alt.y975  = c(rbind(alt.q975,alt.q975))
               mod.y025  = c(rbind(mod.q025,mod.q025))
               mod.y975  = c(rbind(mod.q975,mod.q975))
               keep      = ! (is.na(obs.y025) | is.na(obs.y975))
               iblck     = cumsum(!keep)
               kalt      = ! (is.na(alt.y025) | is.na(alt.y975))
               ibalt     = cumsum(!kalt)
               #----- Split polygons into finite blocks. ----------------------------------#
               x025      = split (x=x025     [keep],f=iblck[keep])
               alt.x025  = split (x=alt.x025 [kalt],f=ibalt[kalt])
               obs.y025  = split (x=obs.y025 [keep],f=iblck[keep])
               alt.y025  = split (x=alt.y025 [kalt],f=ibalt[kalt])
               mod.y025  = split (x=mod.y025 [keep],f=iblck[keep])
               x975      = lapply(X = x025                                   , FUN = rev)
               alt.x975  = lapply(X = alt.x025                               , FUN = rev)
               obs.y975  = lapply(X = split (x=obs.y975 [keep],f=iblck[keep]), FUN = rev)
               alt.y975  = lapply(X = split (x=alt.y975 [kalt],f=ibalt[kalt]), FUN = rev)
               mod.y975  = lapply(X = split (x=mod.y975 [keep],f=iblck[keep]), FUN = rev)

               #----- Plot polygons. ------------------------------------------------------#
               npoly     = length(x025)
               npalt     = length(alt.x025)
               for (y in sequence(npoly)){
                  epolygon( x       = c(    x025[[y]],    x975[[y]])
                          , y       = c(obs.y025[[y]],obs.y975[[y]])
                          , col     = obs.col
                          , angle   =  45
                          , density =  40
                          )#end epolygon
               }#end for
               for (y in sequence(npalt)){
                  epolygon( x       = c(alt.x025[[y]],alt.x975[[y]])
                          , y       = c(alt.y025[[y]],alt.y975[[y]])
                          , col     = alt.col
                          , angle   = -45
                          , density =  40
                          )#end epolygon
               }#end for
               for (y in sequence(npoly)){
                  epolygon( x       = c(    x025[[y]],    x975[[y]])
                          , y       = c(mod.y025[[y]],mod.y975[[y]])
                          , col     = ed22.col
                          , angle   = +90
                          , density =  40
                          )#end epolygon
               }#end poly
               #---------------------------------------------------------------------------#



               #------ Plot the expected values. ------------------------------------------#
               points( x    = fnmean.when
                     , y    = obs.expected
                     , type = "o"
                     , pch  = 16
                     , lwd  = 2
                     , col  = obs.fg
                     )#end points
               points( x    = fnmean.when
                     , y    = alt.expected
                     , type = "o"
                     , pch  = 16
                     , lwd  = 2
                     , col  = alt.fg
                     )#end points
               points( x    = fnmean.when
                     , y    = mod.expected
                     , type = "o"
                     , pch  = 16
                     , lwd  = 2
                     , col  = ed22.fg
                     )#end points
               #---------------------------------------------------------------------------#
            }#end if (plotit)
            #------------------------------------------------------------------------------#

            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main = paste(lyr.title,"\n","N = ",ndat,sep=""),cex.main=0.8*cex.ptsz)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#


         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,ylab=ley,line.ylab=2.5)
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



   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.ftnight)
#------------------------------------------------------------------------------------------#




















#------------------------------------------------------------------------------------------#
#      Plot the time series of fortnightly means.                                          #
#------------------------------------------------------------------------------------------#
if (plot.bp.diel){
   cat(" + Box plot by hour of the day...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      this.compvar  = compvar[[v]]
      this.vnam     = this.compvar$vnam
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = this.compvar$desc
      this.unit     = this.compvar$unit
      this.soilvar  = this.compvar$soilvar
      this.sunvar   = this.compvar$sunvar
      this.col      = this.compvar$col
      this.fg       = this.compvar$fg
      cat("   - ",this.desc,"...","\n")


      simul.ylimit = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      bpdl.list    = list()
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         cat("     * Site :",this.longname,"...","\n")
         #---------------------------------------------------------------------------------#


         #----- Retrieve data. ------------------------------------------------------------#
         obs.now = obser[[this.vnam]]
         nlyr    = ncol(obs.now)
         #---------------------------------------------------------------------------------#



         #----- Find out what to append. --------------------------------------------------#
         if (this.soilvar){
            lyr.key   = obser$slz.key
            lyr.desc  = paste(this.desc,obser$slz.desc,sep="")
            lyr.title = paste(this.longname,"(",toupper(iata),")",obser$slz.desc,sep="")
         }else{
            lyr.key   = ""
            lyr.desc  = this.desc
            lyr.title = paste(this.longname,"(",toupper(iata),")",sep="")
         }#end if
         cnt  = obser$highsun | (! this.sunvar)
         ndat = colSums(is.finite(obs.now[cnt,,drop=FALSE]))
         #---------------------------------------------------------------------------------#


         #------ Plot only if the data are finite. ----------------------------------------#
         obs.list = list()
         mod.list = list()
         for (cc in sequence(nlyr)){

            #------ Split data by lists. --------------------------------------------------#
            obs.list[[cc]]        = split(x=obs.now[,cc],f=obser$hr.idx)
            names(obs.list[[cc]]) = paste("obs",hour.key,sep=".")
            #------------------------------------------------------------------------------#


            #------ Initialise limits for y. ----------------------------------------------#
            ylimit.iata = range(obs.now[,cc],finite=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            mod.list[[cc]]     = list()

            for (s in sequence(nsimul)){
               model   = res[[iata]]$ans[[s]]

               #----- Select and sort the data. -------------------------------------------#
               mod.now = model[[this.vnam]]
               #---------------------------------------------------------------------------#


               #----- Split data into lists. ----------------------------------------------#
               mod.list[[cc]][[s]]        = split(x=mod.now[,cc],f=obser$hr.idx)
               names(mod.list[[cc]][[s]]) = paste("mod",hour.key,sep=".")
               #---------------------------------------------------------------------------#


               #----- Update range. -------------------------------------------------------#
               simul.ylimit[[s]] = range( x      = c( simul.ylimit[[s]]
                                                    , obs.now[,cc]
                                                    , mod.now[,cc]
                                                    )#end c
                                        , finite = TRUE
                                        )#end range
               ylimit.iata       = range( x      = c( ylimit.iata
                                                    , mod.now[,cc]
                                                    )#end c
                                        , finite = TRUE
                                        )#end range
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            names(mod.list[[cc]]) = simul$name
            #------------------------------------------------------------------------------#
         }#end for (cc in sequence(nlyr))
         #---------------------------------------------------------------------------------#


         #----- Save data for this level. -------------------------------------------------#
         bpdl.list[[iata]] = list( obs       = obs.list
                                 , mod       = mod.list
                                 , ylimit    = ylimit.iata
                                 , lyr.key   = lyr.key
                                 , lyr.desc  = lyr.desc
                                 , lyr.title = lyr.title
                                 , ndat      = ndat
                                 , nlyr      = nlyr
                                 , slz       = obser$slz
                                 )#end list
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #      Plot the box plots by scenario.                                               #
      #------------------------------------------------------------------------------------#
      cat("     * Plotting by scenario...","\n")
      for (p in sequence(nsites)){
         iata  = sites$iata[p]
         #----- Get the basic information. ------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         #---------------------------------------------------------------------------------#



         #----- Copy data to local variables. ---------------------------------------------#
         obs.list  = bpdl.list[[iata]]$obs
         mod.list  = bpdl.list[[iata]]$mod
         nlyr      = bpdl.list[[iata]]$nlyr
         ndat      = bpdl.list[[iata]]$ndat
         slz       = bpdl.list[[iata]]$slz
         lyr.key   = bpdl.list[[iata]]$lyr.key
         lyr.desc  = bpdl.list[[iata]]$lyr.desc
         lyr.title = bpdl.list[[iata]]$lyr.title
         ylimit    = pretty.xylim(u=bpdl.list[[iata]]$ylimit,fracexp=0.0,is.log=FALSE)
         lyr.use   = which(ndat > 0)
         #---------------------------------------------------------------------------------#

         #------ Set some common features. ------------------------------------------------#
         letitre = this.longname
         lex     = desc.unit(desc="Time",unit=untab$utc)
         #---------------------------------------------------------------------------------#



         #------ Set some plot defaults. --------------------------------------------------#
         xlimit   = c(0,2*nhour)
         xat      = seq(from=1.5,to=2*nhour-0.5,by=2)
         xlabels  = hour.label
         xgrid    = seq(from=0.5,to=2*nhour+0.5,by=2)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot only those layers with some data.                                      #
         #---------------------------------------------------------------------------------#
         for (cc in lyr.use){

            #------ Set y-axis according to the layer. ------------------------------------#
            ley     = desc.unit(desc=lyr.desc[cc],unit=this.unit)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot the fortnightly means.                                             #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$bp.diel$sites[[this.vnam]]
               fichier = file.path(out.now,paste("bpdiel-",this.vnam,"-",iata
                                                ,lyr.key[cc],".",outform[o],sep="")
                                  )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the device into multiple panels.                                #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(oma=c(1,1,4,0))
               layout(mat = lo.simul$mat)
               #---------------------------------------------------------------------------#



               #----- Loop over all simulations. ------------------------------------------#
               for (s in sequence(nsimul)){

                  #----- Set the lists and colours. ---------------------------------------#
                  bp.list        = c( rbind(obs.list[[cc]],mod.list[[cc]][[s]]))
                  names(bp.list) = c( rbind( names(obs.list[[cc]]     )
                                           , names(mod.list[[cc]][[s]])
                                           )#end rbind
                                    )#end c
                  bp.colour      = rep(c(obs.col,ed22.col),times=nhour)
                  #------------------------------------------------------------------------#



                  #------ Open sub-plot. --------------------------------------------------#
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  if (lo.simul$bottom[s]){
                     axis.rt(side=1,at=xat,labels=xlabels,off=0.05,las=5)
                  }#end if
                  if (lo.simul$left  [s]){
                     axis.rt(side=2,las=1)
                  }#end if
                  abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
                  #------------------------------------------------------------------------#


                  #----- Plot the boxes. --------------------------------------------------#
                  boxplot( x          = bp.list
                         , col        = bp.colour
                         , notch      = TRUE
                         , add        = TRUE
                         , show.names = FALSE
                         , axes       = FALSE
                         )#end boxplot
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s],line=1.0)
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#


               #----- Plot title. ---------------------------------------------------------#
               gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=4.0)
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
         }#end for (cc in lyr.use)
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#






      #------------------------------------------------------------------------------------#
      #      Plot the box plots for default simulation.                                    #
      #------------------------------------------------------------------------------------#
      cat("     * Plotting default simulation...","\n")
      s      = sim.default
      ylimit = pretty.xylim(u=simul.ylimit[[s]],fracexp=0.0,is.log=FALSE)


      #------ Set some common features. ---------------------------------------------------#
      letitre = "Diel distribution"
      lex     = desc.unit(desc="Time",unit=untab$utc)
      ley     = desc.unit(desc=this.desc,unit=this.unit)
      #------------------------------------------------------------------------------------#



      #------ Set some plot defaults. -----------------------------------------------------#
      xlimit   = c(0,2*nhour)
      xat      = seq(from=1.5,to=2*nhour-0.5,by=2)
      xlabels  = hour.label
      xgrid    = seq(from=0.5,to=2*nhour+0.5,by=2)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Plot the fortnightly means.                                                   #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$bp.diel$default
         fichier = file.path(out.now,paste("bpdiel-",this.vnam,"-",simul$name[s],"."
                                          ,outform[o],sep="")
                            )#end file.path
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(1,1,4,0))
         layout(mat = lo.site$mat)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over sites.                                                           #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata  = sites$iata[p]
            #----- Get the basic information. ---------------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #------ Open sub-plot. --------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            if (lo.site$bottom[p]){
               axis.rt(side=1,at=xat,labels=xlabels,off=0.05,las=5)
            }#end if
            if (lo.site$left  [p]){
               axis.rt(side=2,las=1)
            }#end if
            abline(v=xgrid,h=axTicks(2),col=grid.colour,lty="solid")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Plot boxes only if there is anything to plot.                             #
            #------------------------------------------------------------------------------#
            if (iata %in% names(bpdl.list)){
               ndat          = bpdl.list[[iata]]$ndat
               nlyr          = bpdl.list[[iata]]$nlyr
               slz           = bpdl.list[[iata]]$slz
               slz.use       = ifelse(ndat==0,-slz,slz)
               cc            = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
               lyr.key       = bpdl.list[[iata]]$lyr.key  [cc]
               lyr.desc      = bpdl.list[[iata]]$lyr.desc [cc]
               lyr.title     = bpdl.list[[iata]]$lyr.title[cc]
               obs.list      = bpdl.list[[iata]]$obs[[cc]]
               mod.list      = bpdl.list[[iata]]$mod[[cc]][[s]]
               n.list        = sum(sapply(X=obs.list,FUN=is.finite))


               #----- Set the lists and colours. ------------------------------------------#
               bp.list        = c(rbind(obs.list,mod.list))
               names(bp.list) = c(rbind(names(obs.list),names(mod.list)))
               bp.colour      = rep(c(obs.col,ed22.col),times=nhour)
               #---------------------------------------------------------------------------#




               #----- Plot the boxes. -----------------------------------------------------#
               boxplot( x          = bp.list
                      , col        = bp.colour
                      , notch      = TRUE
                      , add        = TRUE
                      , show.names = FALSE
                      , axes       = FALSE
                      )#end boxplot
               #---------------------------------------------------------------------------#
            }#end if (iata %in% names(bpdl.list))


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(lyr.title,"\n","N = ",n.list,sep=""),line=1.0)
            #------------------------------------------------------------------------------#
         }# for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,line.xlab=4.0)
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
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.bp.diel)
#------------------------------------------------------------------------------------------#




















#------------------------------------------------------------------------------------------#
#      Plot the Q-Q plots of daily means.                                                  #
#------------------------------------------------------------------------------------------#
if (plot.qq.dmean){
   cat(" + Find the Q-Q plots of daily means...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      this.compvar  = compvar[[v]]
      this.vnam     = this.compvar$vnam
      this.dmean    = paste("dmean"   ,this.vnam,sep=".")
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = this.compvar$desc
      this.unit     = this.compvar$unit
      this.soilvar  = this.compvar$soilvar
      this.col      = this.compvar$col
      this.fg       = this.compvar$fg
      cat("   - ",this.desc,"...","\n")

      #------ Initialise the global Q-Q plot list. ----------------------------------------#
      qq.list     = list()
      qqn.list    = list()
      qq.limit    = NULL
      qqn.limit   = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      leg.iata    = character(length=nsites); names(leg.iata) = sites$iata
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Loop over sites.                                                              #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         #----- Grab current site code. ---------------------------------------------------#
         iata            = sites$iata[p]
         this.longname   = sites$desc[p]
         #---------------------------------------------------------------------------------#


         #----- Grab data. ----------------------------------------------------------------#
         obser       = eft[[iata]]
         dmean.obser = obser[[this.dmean]]
         ndat        = apply(X = is.finite(dmean.obser), MARGIN=2,FUN=sum)
         nlyr        = length(ndat)
         #---------------------------------------------------------------------------------#



         #----- Normalise the observations. -----------------------------------------------#
         dmean.obser  = split (x=dmean.obser,f=col(dmean.obser))
         mu.obser     = lapply(X = dmean.obser, FUN = mean, na.rm=TRUE)
         sigma.obser  = lapply(X = dmean.obser, FUN = sd  , na.rm=TRUE)
         dnorm.obser  = mapply( FUN      = normalise
                              , x        = dmean.obser
                              , mu       = mu.obser
                              , sigma    = sigma.obser
                              , SIMPLIFY = FALSE
                              )#end mapply
         #---------------------------------------------------------------------------------#



         #----- Find out what to append. --------------------------------------------------#
         if (this.soilvar){
            lyr.key   = obser$slz.key
            lyr.desc  = paste(this.desc,obser$slz.desc,sep="")
            lyr.title = paste(this.longname," (",toupper(iata),")",obser$slz.desc,sep="")
         }else{
            lyr.key   = ""
            lyr.desc  = this.desc
            lyr.title = paste(this.longname," (",toupper(iata),")",sep="")
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Initialise list for this site. --------------------------------------------#
         qq.list [[iata]]  = list( nlyr      = nlyr
                                 , ndat      = ndat
                                 , lyr.key   = lyr.key
                                 , lyr.desc  = lyr.desc
                                 , lyr.title = lyr.title
                                 )#end list
         qqn.list[[iata]]  = qq.list [[iata]]
         slz.use           = ifelse(ndat==0,-obser$slz,obser$slz)
         cc                = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
         leg.iata[p]       = paste(lyr.title[cc]," - N =",ndat[cc],sep="")
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #      Loop over simulations.                                                     #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nsimul)){
            #----- Initialise lists. ------------------------------------------------------#
            qq.list [[iata]][[simul$name[s]]] = list()
            qqn.list[[iata]][[simul$name[s]]] = list()
            #------------------------------------------------------------------------------#



            #----- Select and sort the data. ----------------------------------------------#
            model       = res[[iata]]$ans[[s]]
            dmean.model = split(x=model[[this.dmean]],f=col(model[[this.dmean]]))
            dnorm.model = mapply( FUN      = normalise
                                , x        = dmean.model
                                , mu       = mu.obser
                                , sigma    = sigma.obser
                                , SIMPLIFY = FALSE
                                )#end mapply
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #     Loop over all layers.                                                    #
            #------------------------------------------------------------------------------#
            for (cc in sequence(nlyr)){
               #----- Find the Q-Q plot and copy to the global list. ----------------------#
               qq.now  = qqplot(x=dmean.obser[[cc]],y=dmean.model[[cc]],plot.it=FALSE)
               qqn.now = qqplot(x=dnorm.obser[[cc]],y=dnorm.model[[cc]],plot.it=FALSE)
               qq.list [[iata]][[simul$name[s]]][[cc]] = qq.now
               qqn.list[[iata]][[simul$name[s]]][[cc]] = qqn.now
               #---------------------------------------------------------------------------#



               #----- Update range. -------------------------------------------------------#
               qq.limit       = range(c(qq.limit,qq.now$x,qq.now$y)        ,finite=TRUE)
               qqn.limit[[s]] = range(c(qqn.limit[[s]],qqn.now$x,qqn.now$y),finite=TRUE)
               #---------------------------------------------------------------------------#
            }#end for (cc in sequence(nlyr))
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#





      #------ Set some common features. ---------------------------------------------------#
      letitre = "Q-Q plot of daily means"
      lex     = desc.unit(desc=paste("Observed:" ,this.desc),unit=this.unit)
      ley     = desc.unit(desc=paste("Simulated:",this.desc),unit=this.unit)
      xylimit = pretty.xylim(u=qq.limit,fracexp=0.0,is.log=FALSE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the QQ-plots.                                                            #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$qq.dmean$sites
         fichier = file.path(out.now,paste("qq_dmean-",this.vnam,".",outform[o],sep=""))
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(0,1,4,0))
         layout( mat     = rbind(lo.site$mat.off,rep(1,times=lo.site$ncol))
               , heights = c(rep(6/lo.site$nrow,times=lo.site$nrow),1)
               )#end layout
         #---------------------------------------------------------------------------------#



         #----- Legend. -------------------------------------------------------------------#
         par(mar=c(0.2,0.1,0.1,0.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
         legend ( x       = "bottom"
                , inset   = 0.0
                , legend  = simul$desc
                , fill    = simul$col
                , border  = simul$col
                , ncol    = 3
                , title   = expression(bold("Structure"))
                , cex     = 0.9 * cex.ptsz
                , xpd     = TRUE
                )#end legend
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            #----- Get the basic site information. ----------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #------ Open sub-plot. --------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=xylimit,ylim=xylimit)
            if (lo.site$bottom[p]) axis.rt(side=1,las=5,off=0.06)
            if (lo.site$left  [p]) axis   (side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline(a=0,b=1,col=red.mg,lwd=2,lty="dotdash")
            #------------------------------------------------------------------------------#


            #----- Get Q-Q plot info for this site. ---------------------------------------#
            iata      = sites$iata[p]
            qq.iata   = qq.list[[iata]]
            ndat      = qq.iata$ndat
            nlyr      = qq.iata$nlyr
            cc        = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
            lyr.key   = qq.iata$lyr.key  [cc]
            lyr.desc  = qq.iata$lyr.desc [cc]
            lyr.title = qq.iata$lyr.title[cc]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               #------ Plot the lines and points. -----------------------------------------#
               points( x    = qq.iata[[simul$name[s]]][[cc]]
                     , type = "o"
                     , pch  = 16
                     , cex  = 0.5
                     , col  = simul$col[s]
                     )#end points
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(lyr.title," (n=",ndat[cc],")",sep=""))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #----- Plot title. ---------------------------------------------------------------#
         gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,off.xlab=1/12)
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





      #------------------------------------------------------------------------------------#
      #      Plot the Normalised QQ-plots.                                                 #
      #------------------------------------------------------------------------------------#
      #------ Set some common features. ---------------------------------------------------#
      s       = sim.default
      letitre = "Q-Q plot of daily means (Normalised)"
      lex     = desc.unit(desc=paste("Observed:" ,this.desc),unit=untab$empty)
      ley     = desc.unit(desc=paste("Simulated:",this.desc),unit=untab$empty)
      xylimit = pretty.xylim(u=qqn.limit[[s]],fracexp=0.0,is.log=FALSE)
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$qq.dmean$default
         fichier = file.path(out.now,paste("qqn_dmean-",this.vnam,".",outform[o],sep=""))
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         layout(mat=rbind(2,1),heights=c(5,1))
         #---------------------------------------------------------------------------------#



         #----- Legend. -------------------------------------------------------------------#
         par(mar=c(0.2,3.6,0.1,1.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
         legend ( x       = "bottom"
                , inset   = 0.0
                , legend  = leg.iata
                , fill    = sites$col
                , border  = sites$col
                , ncol    = 2
                , title   = expression(bold("Sites"))
                , cex     = 0.7 * cex.ptsz
                , xpd     = TRUE
                )#end legend
         #---------------------------------------------------------------------------------#



         #------ Open sub-plot. -----------------------------------------------------------#
         par(mar=c(3.6,3.6,3.1,1.1))
         plot.new()
         plot.window(xlim=xylimit,ylim=xylimit)
         axis(side=1)
         axis(side=2,las=1)
         grid(col=grid.colour,lty="dotted")
         abline(a=0,b=1,col=foreground,lwd=2,lty="dotdash")
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            #----- Get the basic site information. ----------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #----- Get Q-Q plot info for this site. ---------------------------------------#
            iata      = sites$iata[p]
            qqn.iata  = qqn.list[[iata]]
            ndat      = qqn.iata$ndat
            nlyr      = qqn.iata$nlyr
            cc        = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
            lyr.key   = qqn.iata$lyr.key  [cc]
            lyr.desc  = qqn.iata$lyr.desc [cc]
            lyr.title = qqn.iata$lyr.title[cc]
            #------------------------------------------------------------------------------#


            #------ Plot the lines and points. --------------------------------------------#
            points( x    = qqn.iata[[simul$name[s]]][[cc]]
                  , type = "o"
                  , pch  = 16
                  , cex  = 0.7
                  , lwd  = 1.5
                  , col  = sites$col[p]
                  )#end points
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #----- Plot title. ---------------------------------------------------------------#
         title(main=letitre,line=1.0)
         title(xlab=lex,ylab=ley,line=2.0)
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
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.qq.dmean)
#------------------------------------------------------------------------------------------#




















#------------------------------------------------------------------------------------------#
#      Plot the density plots of daily means.                                              #
#------------------------------------------------------------------------------------------#
if (plot.density.dmean){
   cat(" + Find the density function of daily means...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      this.compvar  = compvar[[v]]
      this.vnam     = this.compvar$vnam
      this.dmean    = paste("dmean"   ,this.vnam,sep=".")
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = this.compvar$desc
      this.unit     = this.compvar$unit
      this.soilvar  = this.compvar$soilvar
      this.col      = this.compvar$col
      this.fg       = this.compvar$fg
      cat("   - ",this.desc,"...","\n")

      #------ Save all daily means into a list. -------------------------------------------#
      dmean.list = list()
      qlimit     = mapply(FUN=numeric,0*sequence(nsimul),SIMPLIFY=FALSE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      First loop: obtain the daily means.                                           #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         #----- Grab data. ----------------------------------------------------------------#
         iata          = sites$iata[p]
         this.longname = sites$desc[p]
         obser         = eft[[iata]]
         dmean.obser   = obser[[this.dmean]]
         #---------------------------------------------------------------------------------#



         #----- Find out what to append. --------------------------------------------------#
         if (this.soilvar){
            lyr.key   = obser$slz.key
            lyr.desc  = paste(this.desc,obser$slz.desc,sep="")
            lyr.title = paste(this.longname," (",toupper(iata),")",obser$slz.desc,sep="")
         }else{
            lyr.key   = ""
            lyr.desc  = this.desc
            lyr.title = paste(this.longname," (",toupper(iata),")",sep="")
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      First loop over simulations, to obtain the daily means.                    #
         #---------------------------------------------------------------------------------#
         dmean.model = list()
         for (s in sequence(nsimul)){

            #----- Select and grab the data. ----------------------------------------------#
            model            = res[[iata]]$ans[[s]]
            dmean.model[[s]] = split(x=model[[this.dmean]],f=col(model[[this.dmean]]))
            qlimit     [[s]] = range( c(qlimit[[s]],unlist(dmean.obser),dmean.model[[s]]
                                       ,model$soil$soilcp,model$soil$slmsts)
                                    , finite=TRUE
                                    )#end range
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#


         #----- List of daily means. ------------------------------------------------------#
         dmean.list[[iata]] = list( obser     = split(x=dmean.obser,f=col(dmean.obser))
                                  , model     = dmean.model
                                  , ndat      = colSums(is.finite(dmean.obser))
                                  , nlyr      = ncol(dmean.obser)
                                  , slz       = obser$slz
                                  , lyr.key   = lyr.key
                                  , lyr.desc  = lyr.desc
                                  , lyr.title = lyr.title
                                  , soilcp    = model$soil$soilcp [1]
                                  , soilwp    = model$soil$soilwp [1]
                                  , soilfc    = model$soil$sfldcap[1]
                                  , soilpo    = model$soil$slmsts [1]
                                  )#end list
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Standardise the quantiles for the density function.  Always choose a power of  #
      # two for the total length, because it reduces the amount of interpolation.          #
      #------------------------------------------------------------------------------------#
      qlimit.now = pretty.xylim(u=unlist(qlimit),fracexp=0.0,is.log=FALSE)
      qa         = qlimit.now[1]
      qz         = qlimit.now[2]
      quant      = seq(from=qa,to=qz,length.out=n.quant)
      dens.args  = list(n=n.quant,from=qa,to=qz,na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Second loop: obtain the density functions.                                    #
      #------------------------------------------------------------------------------------#
      dlimit    = NULL
      dens.list = NULL
      for (p in sequence(nsites)){
         #----- Initialise list for this site. --------------------------------------------#
         iata    = sites$iata[p]
         dmean   = dmean.list[[iata]]
         ndat    = dmean$ndat
         nlyr    = dmean$nlyr
         slz.use = ifelse(ndat==0,-dmean$slz,dmean$slz)
         cc      = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
         #---------------------------------------------------------------------------------#


         #----- Density function of observations. -----------------------------------------#
         dens.obser = mapply( FUN      = function(x,...) density.safe(x=x,...)$y
                            , x        = dmean$obser
                            , MoreArgs = dens.args
                            , SIMPLIFY = TRUE
                            )#end mapply
         if (any(is.finite(unlist(dmean$obser)))){
            dlimit     = range(c(dlimit,dens.obser[,cc]),finite=TRUE)
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      First loop over simulations, to obtain the daily means.                    #
         #---------------------------------------------------------------------------------#
         dens.model = list()
         for (s in sequence(nsimul)){
            #----- Find daily mean (but only for those days with full record). ------------#
            dens.model[[s]] = mapply( FUN      = function(x,...) density.safe(x=x,...)$y
                                    , x        = dmean$model[[s]]
                                    , MoreArgs = dens.args
                                    , SIMPLIFY = TRUE
                                    )#end mapply
            if (any(is.finite(unlist(dmean$model[[s]])))){
               dlimit          = range(c(dlimit,dens.model[[s]][,cc]),finite=TRUE)
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#


         #----- List of daily means. ------------------------------------------------------#
         dens.list[[iata]] = list( obser     = dens.obser
                                 , model     = dens.model
                                 , ndat      = dmean$ndat
                                 , nlyr      = dmean$nlyr
                                 , slz       = dmean$slz
                                 , lyr.key   = dmean$lyr.key
                                 , lyr.desc  = dmean$lyr.desc
                                 , lyr.title = dmean$lyr.title
                                 , soilwp    = dmean$soilwp
                                 , soilfc    = dmean$soilfc
                                 , soilpo    = dmean$soilpo
                                 )#end list
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      dlimit = pretty.xylim(u=dlimit,fracexp=0.0,is.log=FALSE)
      #------------------------------------------------------------------------------------#








      #------ Set some common features. ---------------------------------------------------#
      letitre = "Density function: daily means"
      lex     = desc.unit(desc=this.desc,unit=this.unit)
      ley     = desc.unit(desc="Density function",unit=untab$empty)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the density function plots.                                              #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$density.dmean$sites
         fichier = file.path(out.now,paste("density-",this.vnam,".",outform[o],sep=""))
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         if (density.legend){
            #----- Include legend beneath the plots. --------------------------------------#
            par(par.user)
            par(oma=c(0,1,4,0))
            layout( mat     = rbind(lo.site$mat.off,rep(1,times=lo.site$ncol))
                  , heights = c(rep(6/lo.site$nrow,times=lo.site$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#




            #----- Legend. ----------------------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , fill    = simul$col
                   , border  = simul$col
                   , ncol    = 3
                   , title   = expression(bold("Structure"))
                   , cex     = 0.9 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#

         }else{
            #----- Don't add legend beneath the plots. ------------------------------------#
            par(par.user)
            par(oma=c(0,1,4,0))
            layout(mat=lo.site$mat)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            #----- Get the basic site information. ----------------------------------------#
            iata          = sites$iata[p]
            this.longname = sites$desc[p]
            #------------------------------------------------------------------------------#



            #------ Open sub-plot. --------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=qlimit.now,ylim=dlimit)
            if (lo.site$bottom[p]) axis(side=1,las=1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            #------------------------------------------------------------------------------#


            #----- Get density plot info for this site. -----------------------------------#
            iata      = sites$iata[p]
            dens.iata = dens.list[[iata]]
            ndat      = dens.iata$ndat
            nlyr      = dens.iata$nlyr
            slz.use   = ifelse(ndat==0,-dens.iata$slz,dens.iata$slz)
            cc        = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
            lyr.key   = dens.iata$lyr.key  [cc]
            lyr.desc  = dens.iata$lyr.desc [cc]
            lyr.title = dens.iata$lyr.title[cc]
            #------------------------------------------------------------------------------#



            #----- If this is soil moisture, add the lines with special thresholds. -------#
            if (this.vnam %in% c("soil.water")){
               abline(v = dens.iata$soilcp, col = "#990F0F",lwd=2.0,lty="dotdash")
               abline(v = dens.iata$soilwp, col = "#E65C17",lwd=2.0,lty="dotdash")
               abline(v = dens.iata$soilfc, col = "#2996CC",lwd=2.0,lty="dotdash")
               abline(v = dens.iata$soilpo, col = "#3B24B3",lwd=2.0,lty="dotdash")
            }#end if
            #------------------------------------------------------------------------------#




            #----- Plot the density function of observations. -----------------------------#
            lines ( x    = quant
                  , y    = dens.iata$obser[,cc]
                  , type = "l"
                  , lwd  = 2.5
                  , col  = foreground
                  )#end points
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               #------ Plot the lines. ----------------------------------------------------#
               lines ( x    = quant
                     , y    = dens.iata$model[[s]][,cc]
                     , type = "l"
                     , lwd  = 2.0
                     , col  = simul$col[s]
                     )#end points
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(lyr.title,"\n","N = ",ndat[cc],sep=""),line=0.25)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #----- Plot title. ---------------------------------------------------------------#
         if (density.legend){
            gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,off.xlab=1/12)
         }else{
            gtitle(main=letitre,xlab=lex,ylab=ley,line.xlab=4.2,line.ylab=2.5)
         }#end if
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






      #------------------------------------------------------------------------------------#
      #     Standardise the quantiles for the density function.  Always choose a power of  #
      # two for the total length, because it reduces the amount of interpolation.          #
      #------------------------------------------------------------------------------------#
      s          = sim.default
      qlimit.now = pretty.xylim(u=unlist(qlimit[[s]]),fracexp=0.0,is.log=FALSE)
      qa    = qlimit.now[1]
      qz    = qlimit.now[2]
      quant = seq(from=qa,to=qz,length.out=n.quant)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Second loop: obtain the density functions.                                    #
      #------------------------------------------------------------------------------------#
      dlimit    = NULL
      dens.list = NULL
      for (p in sequence(nsites)){
         #----- Initialise list for this site. --------------------------------------------#
         iata    = sites$iata[p]
         dmean   = dmean.list[[iata]]
         ndat    = dmean$ndat
         nlyr    = dmean$nlyr
         slz.use = ifelse(ndat==0,-dmean$slz,dmean$slz)
         cc      = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
         #---------------------------------------------------------------------------------#





         #----- Density function of observations. -----------------------------------------#
         dens.obser = mapply( FUN      = function(x,...) density.safe(x=x,...)$y
                            , x        = dmean$obser
                            , MoreArgs = dens.args
                            , SIMPLIFY = TRUE
                            )#end mapply
         if (any(is.finite(unlist(dmean$obser)))){
            dlimit     = range(c(dlimit,dens.obser[,cc]),finite=TRUE)
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      First loop over simulations, to obtain the daily means.                    #
         #---------------------------------------------------------------------------------#
         dens.model = mapply( FUN      = function(x,...) density.safe(x=x,...)$y
                            , x        = dmean$model[[s]]
                            , MoreArgs = dens.args
                            , SIMPLIFY = TRUE
                            )#end mapply
         if (any(is.finite(unlist(dmean$model[[s]])))){
            dlimit     = range(c(dlimit,dens.model[,cc]),finite=TRUE)
         }#end if
         #---------------------------------------------------------------------------------#


         #----- List of daily means. ------------------------------------------------------#
         dens.list[[iata]] = list( obser     = dens.obser
                                 , model     = dens.model
                                 , ndat      = dmean$ndat
                                 , nlyr      = dmean$nlyr
                                 , slz       = dmean$slz
                                 , lyr.key   = dmean$lyr.key
                                 , lyr.desc  = dmean$lyr.desc
                                 , lyr.title = dmean$lyr.title
                                 , soilcp    = dmean$soilcp
                                 , soilwp    = dmean$soilwp
                                 , soilfc    = dmean$soilfc
                                 , soilpo    = dmean$soilpo
                                 )#end list
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      dlimit = pretty.xylim(u=dlimit,fracexp=0.0,is.log=FALSE)
      #------------------------------------------------------------------------------------#








      #------ Set some common features. ---------------------------------------------------#
      letitre = "Density function: daily means"
      lex     = desc.unit(desc=this.desc,unit=this.unit)
      ley     = desc.unit(desc="Density function",unit=untab$empty)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the density function plots.                                              #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$density.dmean$default
         fichier = file.path(out.now,paste("density-",this.vnam,".",outform[o],sep=""))
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
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Split the device into multiple panels.                                      #
         #---------------------------------------------------------------------------------#
         par(par.user)
         par(oma=c(0,1,4,0))
         layout(mat=lo.site$mat)
         #---------------------------------------------------------------------------------#



         #----- Loop over all sites. ------------------------------------------------------#
         for (p in sequence(nsites)){
            #----- Get density plot info for this site. -----------------------------------#
            iata      = sites$iata[p]
            dens.iata = dens.list[[iata]]
            ndat      = dens.iata$ndat
            nlyr      = dens.iata$nlyr
            slz.use   = ifelse(ndat==0,-dens.iata$slz,dens.iata$slz)
            cc        = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
            lyr.key   = dens.iata$lyr.key  [cc]
            lyr.desc  = dens.iata$lyr.desc [cc]
            lyr.title = dens.iata$lyr.title[cc]
            #------------------------------------------------------------------------------#


            #------ Open sub-plot. --------------------------------------------------------#
            par(mar=lo.site$mar[p,])
            plot.new()
            plot.window(xlim=qlimit.now,ylim=dlimit)
            if (lo.site$bottom[p]) axis(side=1,las=1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            #------------------------------------------------------------------------------#



            #----- If this is soil moisture, add the lines with special thresholds. -------#
            if (this.vnam %in% "soil.water"){
               abline(v = dens.iata$soilcp, col = "#990F0F",lwd=2.0,lty="dotdash")
               abline(v = dens.iata$soilwp, col = "#E65C17",lwd=2.0,lty="dotdash")
               abline(v = dens.iata$soilfc, col = "#2996CC",lwd=2.0,lty="dotdash")
               abline(v = dens.iata$soilpo, col = "#3B24B3",lwd=2.0,lty="dotdash")
            }#end if
            #------------------------------------------------------------------------------#


            #----- Plot the density function of observations and predictions. -------------#
            lines ( x    = quant
                  , y    = dens.iata$obser[,cc]
                  , type = "l"
                  , lwd  = 2.5
                  , col  = obs.fg
                  )#end points
            lines ( x    = quant
                  , y    = dens.iata$model[,cc]
                  , type = "l"
                  , lwd  = 2.5
                  , col  = ed22.col
                  )#end points
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(lyr.title,"\n","N = ",ndat[cc],sep=""),line=0.25)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#



         #----- Plot title. ---------------------------------------------------------------#
         if (density.legend){
            gtitle(main=letitre,xlab=lex,ylab=ley,line.ylab=2.5,off.xlab=1/12)
         }else{
            gtitle(main=letitre,xlab=lex,ylab=ley,line.xlab=4.2,line.ylab=2.5)
         }#end if
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


   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.density.dmean)
#------------------------------------------------------------------------------------------#


























#------------------------------------------------------------------------------------------#
#      Plot the spider web for all sites and all variables, for statistics that may can be #
# plotted in a spider web (i.e., positive defined).                                        #
#------------------------------------------------------------------------------------------#
if (plot.spider){
   cat(" + Plot the spider web diagrams for all sites and all variables...","\n")
   good.loop = which(unlist(sapply(X=good,FUN=c)["spider",]))
   for (g in good.loop){
      #---- Copy structured variables to convenient scratch scalars. ----------------------#
      this.good = good[[g]]$vnam
      desc.good = good[[g]]$desc
      norm.good = good[[g]]$normalise
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Load all data to an array.                                                    #
      #------------------------------------------------------------------------------------#
      cat("   - ",desc.good,"...","\n")
      web = array( dim      = c   (nsimul,nsites,ncompvar,ndiel,nseason)
                 , dimnames = list(simul.key,sites$iata,compvar.key,diel.key,season.key)
                 )#end array 
      for (v in sequence(ncompvar)){
         this.vnam     = compvar[[v]]$vnam
         this.soilvar  = compvar[[v]]$soilvar



         #---------------------------------------------------------------------------------#
         #     Loop over all sites.                                                        #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata         = sites$iata[p]
            longname     = sites$desc[p]
            obser        = eft[[iata]]
            obs.now      = obser[[this.vnam]]
            nlyr         = ncol(obs.now)
            ndat         = colSums(is.finite(obs.now))
            slz.use      = ifelse(ndat==0,-obser$slz,obser$slz)
            cc           = pmin(nlyr,which.min(abs(slz.use-slz.reference)))


            #------------------------------------------------------------------------------#
            #     Grab the data for this simulation.                                       #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               comp.now     = res[[iata]]$sim[[s]][[this.vnam]]
               good.now     = comp.now[[this.good]][cc,,]
               sfac.now     = 1. + norm.good * ( sqrt(comp.now$obs.moment[cc,,,2]) - 1 )
               sfac.now     = ifelse(sfac.now == 0.,NA,1/sfac.now)
               web[s,p,v,,] = abs(good.now) * sfac.now
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end for (p in 1:nsites)
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:ncompvar)
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Plot the spider webs by diel and variable.                                     #
      #------------------------------------------------------------------------------------#
      for (d in sequence(ndiel)){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Webs by site (all variables).                                               #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata     = sites$iata[p]

            letitre = paste(desc.good," - ",sites$desc[p],"\n",diel.desc[d],sep="")

            if (any(is.finite(web[,p,,d,nseason]))){
               v.sel = is.finite(colSums(web[,p,,d,nseason]))

               if (this.good %in% "sw.stat"){
                  web.range = c(0,1)
               }else{
                  web.range = range(c(0,web[,p,v.sel,d,nseason]),na.rm=TRUE)
               }#end if
               if (ptsz <= 11){
                  web.lim   = pretty(web.range,n=5)
               }else if (ptsz <= 14){
                  web.lim   = pretty(web.range,n=4)
               }else{
                  web.lim   = pretty(web.range,n=3)
               }#end if

               #---------------------------------------------------------------------------#
               #     Webs by variable (all sites).                                         #
               #---------------------------------------------------------------------------#
               for (o in sequence(nout)){
                  #----- Make the file name. ----------------------------------------------#
                  out.web = out[[outform[o]]]$spider[[diel.key[d]]]$sites
                  fichier   = file.path(out.web,paste("spider-",this.good,"-",iata
                                               ,"-",diel.key[d],".",outform[o],sep="")
                                       )#end file.path
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
                  #     Split the window into 3, and add site and simulation legends at    #
                  # the bottom.                                                            #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(oma = c(0.2,0,3.0,0))
                  layout(mat = rbind(2,1),height = c(18,4))
                  #------------------------------------------------------------------------#




                  #----- Legend: the simulations. -----------------------------------------#
                  par(mar=c(0.1,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = simul$desc
                         , fill    = simul$col
                         , border  = simul$col
                         , ncol    = 2
                         , title   = expression(bold("Structure"))
                         , cex     = 0.75 * cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the spider web.                                               #
                  #------------------------------------------------------------------------#
                  radial.flex( lengths          = web[,p,v.sel,d,nseason]
                             , labels           = as.expression(compvar.sym[v.sel])
                             , lab.col          = foreground
                             , lab.bg           = background
                             , radlab           = FALSE
                             , start            = 90
                             , clockwise        = TRUE
                             , rp.type          = "p"
                             , label.prop       = 1.15 * max(1,sqrt(ptsz / 14))
                             , main             = ""
                             , line.col         = simul$col
                             , lty              = simul$lty
                             , lwd              = 3.0
                             , show.grid        = TRUE
                             , show.grid.labels = 4
                             , show.radial.grid = TRUE
                             , grid.col         = grid.colour
                             , radial.lim       = web.lim
                             , radial.col       = foreground
                             , radial.bg        = background
                             , poly.col         = NA
                             , mar              = c(2,1,2,1)+0.1
                             , cex.lab          = 0.5
                             )#end radial.plot
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  par(las=0)
                  mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
                  #------------------------------------------------------------------------#


                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  dummy = clean.tmp()
                  #------------------------------------------------------------------------#

               }#end for (o in 1:nout)
               #---------------------------------------------------------------------------#
            }#end if (any(is.finite(web[,,v,d,nseason])))
            #------------------------------------------------------------------------------#
         }#end for (v in 1:ncompvar)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Webs by variable (all sites).                                               #
         #---------------------------------------------------------------------------------#
         for (v in sequence(ncompvar)){
            this.vnam     = compvar[[v]]$vnam
            this.desc     = compvar[[v]]$desc

            letitre = paste(desc.good," - ",this.desc,"\n",diel.desc[d],sep="")

            if (any(is.finite(web[,,v,d,nseason]))){
               p.sel = is.finite(colSums(web[,,v,d,nseason]))

               if (this.good %in% "sw.stat"){
                  web.range = c(0,1)
               }else{
                  web.range = range(c(0,web[,p.sel,v,d,nseason]),na.rm=TRUE)
               }#end if
               web.lim   = pretty(web.range,n=4)

               #---------------------------------------------------------------------------#
               #     Webs by variable (all sites).                                         #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Make the file name. ----------------------------------------------#
                  out.web = out[[outform[o]]]$spider[[diel.key[d]]]$variables
                  fichier   = file.path(out.web,paste("spider-",this.good,"-",this.vnam,"-"
                                                     ,diel.key[d],".",outform[o],sep=""))
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
                  #     Split the window into 3, and add site and simulation legends at    #
                  # the bottom.                                                            #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(oma = c(0.2,0,3.0,0))
                  layout(mat = rbind(2,1),height = c(6.0,1.0))
                  #------------------------------------------------------------------------#




                  #----- Legend: the simulations. -----------------------------------------#
                  par(mar=c(0.1,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = simul$desc
                         , fill    = simul$col
                         , border  = simul$col
                         , ncol    = 3
                         , title   = expression(bold("Structure"))
                         , pt.cex  = simul$cex
                         , cex     = 0.75 * cex.ptsz
                         )#end legend
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the spider web.                                               #
                  #------------------------------------------------------------------------#
                  radial.flex( lengths          = web[,p.sel,v,d,nseason]
                             , labels           = toupper(sites$iata[p.sel])
                             , radlab           = FALSE
                             , start            = 90
                             , clockwise        = TRUE
                             , rp.type          = "p"
                             , main             = ""
                             , line.col         = simul$col
                             , lty              = simul$lty
                             , lwd              = 3.0
                             , show.grid        = TRUE
                             , show.grid.labels = 4
                             , show.radial.grid = TRUE
                             , grid.col         = grid.colour
                             , radial.lim       = web.lim
                             , poly.col         = NA
                             , mar              = c(2,1,2,1)+0.1
                             , cex.lab          = 0.5
                             )#end radial.plot
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  par(las=0)
                  mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
                  #------------------------------------------------------------------------#


                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  dummy = clean.tmp()
                  #------------------------------------------------------------------------#

               }#end for (o in 1:nout)
               #---------------------------------------------------------------------------#
            }#end if (any(is.finite(web[,,v,d,nseason])))
            #------------------------------------------------------------------------------#
         }#end for (v in 1:ncompvar)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end for (d in 1:ndiel)
      #------------------------------------------------------------------------------------#




      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Webs by variable (default simulation -- diel).                                 #
      #------------------------------------------------------------------------------------#
      s     = sim.default
      d.sel = diel.default
      for (v in sequence(ncompvar)){
         this.vnam     = compvar[[v]]$vnam
         this.desc     = compvar[[v]]$desc

         letitre = paste(desc.good," - ",this.desc,sep="")

         if (any(is.finite(web[s,,v,d.sel,nseason]))){
            #----- Select only the sites with all diel statistics. ------------------------#
            p.sel  = rowSums(is.finite(web[s,,v,d.sel,nseason])) > 0
            dcount = colSums(is.finite(web[s,,v,,nseason]))
            d.sel  = diel.default  & dcount == max(dcount)
            #------------------------------------------------------------------------------#

            if (this.good %in% "sw.stat"){
               web.range = c(0,1)
            }else{
               web.range = range(c(0,web[s,p.sel,v,d.sel,nseason]),finite=TRUE)
            }#end if
            web.lim   = pretty(web.range,n=4)

            #------------------------------------------------------------------------------#
            #     Webs by variable (all sites).                                            #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Make the file name. -------------------------------------------------#
               out.web = out[[outform[o]]]$spider$default.var
               fichier   = file.path(out.web,paste("spider-",this.good,"-",this.vnam,"-"
                                                  ,simul$name[s],".",outform[o],sep=""))
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into 3, and add site and simulation legends at       #
               # the bottom.                                                               #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,0,3.0,0))
               layout(mat = rbind(2,1),height = c(6.0,1.0))
               #---------------------------------------------------------------------------#




               #----- Legend: the simulations. --------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = diel.desc[d.sel]
                      , fill    = diel.col [d.sel]
                      , border  = diel.col [d.sel]
                      , ncol    = 3
                      , pt.cex  = simul$cex
                      , cex     = 0.75 * cex.ptsz
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the spider web.                                                  #
               #---------------------------------------------------------------------------#
               web.now = web[s,,v,,nseason]
               web.now = t(web.now[p.sel,d.sel,drop=FALSE])
               radial.flex( lengths          = web.now
                          , labels           = toupper(sites$iata[p.sel])
                          , radlab           = FALSE
                          , start            = 90
                          , clockwise        = TRUE
                          , rp.type          = "p"
                          , main             = ""
                          , line.col         = diel.col[d.sel]
                          , lty              = "solid"
                          , lwd              = 3.0
                          , show.grid        = TRUE
                          , show.grid.labels = 4
                          , show.radial.grid = TRUE
                          , grid.col         = grid.colour
                          , radial.lim       = web.lim
                          , poly.col         = NA
                          , mar              = c(2,1,2,1)+0.1
                          , cex.lab          = 0.5
                          )#end radial.plot
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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

            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end if (any(is.finite(web[,,v,d,nseason])))
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:ncompvar)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#







      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Webs by diel (all variables and sites).                                        #
      #------------------------------------------------------------------------------------#
      s         = sim.default
      ee        = nseason
      v.inc     = c("gpp","nep","reco","rshortup","rlongup","ustar","hflxca","wflxca")
      v.use     = compvar.key %in% v.inc
      diel.loop = sequence(ndiel)
      d.sel     = which(diel.key %in% c("all.hrs","dmean"))
      for (d in sequence(ndiel)){
         letitre = paste(desc.good," - ",diel.desc[d],sep="")


         #---------------------------------------------------------------------------------#
         #     Check whether there is anything to plot.                                    #
         #---------------------------------------------------------------------------------#
         if (any(is.finite(web[s,,v.use,d,nseason]))){

            v.keep  = which(v.use & is.finite(colSums(web[s,,,d,ee])))
            v.order = match(v.inc,compvar.key[v.keep])
            v.order = v.order[! is.na(v.order)]
            v.sel   = v.keep[v.order]


            if (this.good %in% "sw.stat"){
               web.range = c(0,1)
            }else{
               if (d %in% d.sel){
                  web.range = range(c(0,web[s,,v.sel,d,ee]),finite=TRUE)
               }else{
                  web.range = range(c(0,web[s,,v.sel,d,ee]),finite=TRUE)
               }#end if
            }#end if
            if (ptsz <= 11){
               web.lim   = pretty(web.range,n=5)
            }else if (ptsz <= 14){
               web.lim   = pretty(web.range,n=4)
            }else{
               web.lim   = pretty(web.range,n=3)
            }#end if


            #------------------------------------------------------------------------------#
            #     Format loop.  We now plot all sites and variables.                       #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.web = out[[outform[o]]]$spider$default.site
               fichier   = file.path(out.web,paste("spider-",this.good,"-",diel.key[d]
                                                  ,"-",simul.key[s],".",outform[o],sep="")
                                    )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into 3, and add site and simulation legends at       #
               # the bottom.                                                               #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,0,3.0,0))
               layout(mat = rbind(2,1),height = c(18,4))
               #---------------------------------------------------------------------------#




               #----- Legend: the simulations. --------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = paste(sites$desc," (",toupper(sites$iata),")",sep="")
                      , col     = sites$col
                      , lwd     = 3.0
                      , ncol    = 3
                      , title   = expression(bold("Sites"))
                      , cex     = 0.75 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Plot the spider webs for all variables.                               #
               #---------------------------------------------------------------------------#
               radial.flex( lengths          = web[s,,v.sel,d,ee]
                          , labels           = as.expression(compvar.sym[v.sel])
                          , lab.col          = foreground
                          , lab.bg           = background
                          , radlab           = FALSE
                          , start            = 90
                          , clockwise        = TRUE
                          , rp.type          = "p"
                          , label.prop       = 1.15 * max(1,sqrt(ptsz / 14))
                          , main             = ""
                          , line.col         = sites$col
                          , lty              = "solid"
                          , lwd              = 3.0
                          , show.grid        = TRUE
                          , show.grid.labels = 2
                          , show.radial.grid = TRUE
                          , grid.col         = grid.colour
                          , radial.lim       = web.lim
                          , radial.col       = foreground
                          , radial.bg        = background
                          , poly.col         = NA
                          , mar              = c(2,1,2,1)+0.1
                          , cex.lab          = 0.5
                          )#end radial.plot
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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

            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end if (any(is.finite(web[s,,,d,nseason])))
         #---------------------------------------------------------------------------------#
      }#end for (d in sequence(ndiel))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   }#end for (g in good.loop)
   #---------------------------------------------------------------------------------------#
}#end if (plot.spider)
#------------------------------------------------------------------------------------------#


















#------------------------------------------------------------------------------------------#
#         Plot the Skill and Taylor diagrams.                                              #
#------------------------------------------------------------------------------------------#
if (plot.skill.taylor){
   cat (" + Plot cross-model and cross-site diagrams (Skill and Taylor)...","\n")
   for (v in sequence(ncompvar)){
      #----- Copy the variable information. -----------------------------------------------#
      this.vnam     = compvar[[v]]$vnam
      this.dmean    = paste("dmean"   ,this.vnam,sep=".")
      this.fnmean   = paste("fnmean"  ,this.vnam,sep=".")
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = compvar[[v]]$desc
      this.unit     = compvar[[v]]$unit
      this.sun      = compvar[[v]]$sunvar
      this.soilvar  = compvar[[v]]$soilvar
      cat("   - ",this.desc,"...","\n")
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Loop over all parts of the day.                                               #
      #------------------------------------------------------------------------------------#
      for (d in sequence(ndiel)){
         cat("     * ",diel.desc[d],"...","\n")


         #---------------------------------------------------------------------------------#
         #      Loop over all sites, normalise data and create the vector for the model.   #
         #---------------------------------------------------------------------------------#
         obs.diel    = list()
         mod.diel    = list()
         cnt.diel    = rep(x=NA,times=nsites); names(cnt.diel) = sites$iata
         bias.range  = NULL
         sigma.range = NULL
         for (p in sequence(nsites)){
            #----- Grab data. -------------------------------------------------------------#
            iata  = sites$iata[p]
            obs   = eft[[iata]]
            nwhen = obs$nwhen
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Decide which data set to use.                                           #
            #------------------------------------------------------------------------------#
            ndat    = colSums(is.finite(obs[[this.vnam]]))
            nlyr    = length(ndat)
            slz.use = ifelse(ndat==0,-obs$slz,obs$slz)
            cc      = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
            if (d %in% diel.fnmean){
               this.obs = obs[[this.fnmean]][,cc]
               sel      = is.finite(this.obs)
               this.obs = this.obs[sel]
            }else if (d %in% diel.dmean){
               this.obs = obs[[this.dmean]][,cc]
               sel      = is.finite(this.obs)
               this.obs = this.obs[sel]
            }else{
               #----- Select this diel (or everything for all day). -----------------------#
               d.sel    = (obs$diel == d | d %in% diel.all.hrs)
               s.sel    = obs$highsun | (! this.sun)
               o.sel    = is.finite(obs[[this.vnam]][,cc])
               sel      = d.sel & s.sel & o.sel
               sel      = ifelse(is.na(sel),FALSE,sel)
               this.obs = obs[[this.vnam]][sel,cc]
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Find the standard deviation of this observation.  Skip the site if       #
            # everything is zero.                                                          #
            #------------------------------------------------------------------------------#
            comp         = res[[iata]]$sim[[simul.key[1]]][[this.vnam]]
            sdev.obs.now = sqrt(comp$obs.moment[cc,d,nseason,2])
            sel          = sel & sdev.obs.now %>% 0
            n.sel        = sum(sel)
            #------------------------------------------------------------------------------#



            #----- Copy the observed data. ------------------------------------------------#
            obs.diel[[iata]] = this.obs
            mod.diel[[iata]] = matrix(ncol=nsimul,nrow=n.sel,dimnames=list(NULL,simul.key))
            #------------------------------------------------------------------------------#



            #----- Copy the modelled data, and update ranges. -----------------------------#
            if (any(sel)){
               for (s in sequence(nsimul)){


                  #------------------------------------------------------------------------#
                  #      Decide which variable to use depending on the diel.               #
                  #------------------------------------------------------------------------#
                  mod  = res[[iata]]$ans[[simul.key[s]]]
                  if (d %in% diel.fnmean){
                     this.mod = mod[[this.fnmean]][sel,cc]
                  }else if (d %in% diel.dmean){
                     this.mod = mod[[this.dmean]][sel,cc]
                  }else{
                     this.mod = mod[[this.vnam]][sel,cc]
                  }#end if
                  mod.diel[[iata]][,s] = this.mod
                  #------------------------------------------------------------------------#


                  #----- Check number of valid entries. -----------------------------------#
                  if (! is.null(cnt.diel[[iata]])){
                     cnt.diel[[iata]] = sum(is.finite(this.obs))
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Find the normalised bias and model standard deviation. -----------#
                  comp         = res[[iata]]$sim[[simul.key[s]]][[this.vnam]]
                  sdev.obs.now = sqrt(comp$obs.moment[cc,d,nseason,2])
                  bias.now     = comp$bias [cc,d,nseason] / sdev.obs.now
                  sigma.now    = comp$sigma[cc,d,nseason] / sdev.obs.now
                  bias.range   = c(bias.range ,bias.now   )
                  sigma.range  = c(sigma.range,sigma.now  )
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#
            }#end if (any(sel))
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Plot Taylor and skill plots only if there is anything to plot.              #
         #---------------------------------------------------------------------------------#
         ok.taylor.skill = (  length(unlist(obs.diel)) > 0  && any(cnt.diel > 0)           
                           && any(is.finite(bias.range))    && any(is.finite(sigma.range)))
         if (ok.taylor.skill){
            #---- Fix ranges. -------------------------------------------------------------#
            xy.range    = 1.04 * max(abs(c(bias.range,sigma.range)),na.rm=TRUE)
            bias.range  = 1.04 * xy.range  * c(-1,1)
            sigma.range = 1.04 * xy.range  * c( 1,0)
            r2.range    = range(1-xy.range^2,1)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Calculate the size for the points in the Skill and Taylor diagrams.      #
            # Make it proportional to the number of points used to evaluate each place.    #
            #------------------------------------------------------------------------------#
            st.cnt.min = min (cnt.diel[cnt.diel > 0] , na.rm = TRUE)
            st.cnt.max = max (cnt.diel[cnt.diel > 0] , na.rm = TRUE)
            st.cnt.med = round(mean(c(st.cnt.min,st.cnt.max)))
            cex.diel   = pmax( st.cex.min, ( st.cex.min + ( st.cex.max  - st.cex.min )
                                                        * (    cnt.diel - st.cnt.min )
                                                        / ( st.cnt.max  - st.cnt.min )
                                           )#end cex.diel
                             )#end pmax
            lwd.diel   = pmax( st.lwd.min, ( st.lwd.min + ( st.lwd.max  - st.lwd.min )
                                                        * (    cnt.diel - st.cnt.min )
                                                        / ( st.cnt.max  - st.cnt.min )
                                           )#end cex.diel
                            )#end pmax
            st.cex.med = mean(c(st.cex.min,st.cex.max))
            #------------------------------------------------------------------------------#


            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #     Skill plot.                                                              #
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Plot title.                                                              #
            #------------------------------------------------------------------------------#
            letitre = paste(" Skill diagram - ",this.desc,"\n",diel.desc[d],sep="")
            cat("       - Skill","\n")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all formats.                                                  #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.skill = out[[outform[o]]]$skill[[diel.key[d]]]
               fichier   = file.path(out.skill,paste("skill-",this.vnam,"-",diel.key[d]
                                                    ,".",outform[o],sep="")
                                    )#end file.path
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into 3, and add site and simulation legends at the   #
               # bottom.                                                                   #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,3.0,0))
               layout(mat = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3)),height = c(5.0,1.0))
               #---------------------------------------------------------------------------#




               #----- Legend: the sites. --------------------------------------------------#
               par(mar=c(0.2,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = toupper(sites$iata)
                      , col     = foreground
                      , pt.bg   = foreground
                      , pch     = sites$pch
                      , ncol    = min(4,pretty.box(nsites)$ncol)
                      , title   = expression(bold("Sites"))
                      , pt.cex  = st.cex.med
                      , cex     = 1.1 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#




               #----- Legend: the counts. -------------------------------------------------#
               par(mar=c(0.2,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = c(st.cnt.min,st.cnt.med,st.cnt.max)
                      , col     = foreground
                      , pt.bg   = foreground
                      , pch     = 15
                      , ncol    = 1
                      , title   = expression(bold("Number Obs."))
                      , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                      , cex     = 1.0 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#




               #----- Legend: the simulations. --------------------------------------------#
               par(mar=c(0.2,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = simul$desc
                      , fill    = simul$col
                      , border  = simul$col
                      , ncol    = 2
                      , title   = expression(bold("Structure"))
                      , cex     = 1.0 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Loop over sites.                                                      #
               #---------------------------------------------------------------------------#
               myskill = NULL
               for (p in sequence(nsites)){
                  iata = sites$iata[p]

                  #----- Skip the site if there is no data. -------------------------------#
                  ok.iata = (  length(obs.diel[[iata]]) > 0
                            && any(is.finite(obs.diel[[iata]])) )
                  ok.iata = ok.iata && ( ! is.na(ok.iata))
                  if (ok.iata){
                     #---------------------------------------------------------------------#
                     #     Initialise or update the skill plot.                            #
                     #---------------------------------------------------------------------#
                     myskill = skill.plot( obs           = obs.diel[[iata]]
                                         , obs.options   = list( col = foreground
                                                               , cex = 2.0
                                                               )#end list
                                         , mod           = mod.diel[[iata]]
                                         , mod.options   = list( col = simul$col
                                                               , bg  = simul$col
                                                               , pch = sites$pch[p]
                                                               , cex = cex.diel [p]
                                                               , lty = "solid"
                                                               , lwd = lwd.diel [p]
                                                               )#end list
                                         , main           = ""
                                         , bias.lim       = bias.range
                                         , r2.lim         = r2.range
                                         , r2.options     = list( col = grid.colour)
                                         , nobias.options = list( col = khaki.mg   )
                                         , rmse.options   = list( col = orange.mg
                                                                , lty = "dotdash"
                                                                , lwd = 1.2
                                                                , bg  = background
                                                                )#end list
                                         , cex.xyzlab     = 1.4
                                         , cex.xyzat      = 1.4
                                         , skill          = myskill
                                         , normalise      = TRUE
                                         , mar            = c(5,4,4,3)+0.1
                                         )#end skill.plot
                     #---------------------------------------------------------------------#
                  }#end if (length(obs.diel[[iata]] > 0)
                  #------------------------------------------------------------------------#
               }#end for (p in 1:nsites)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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
            }#end for (o in 1:nout) 
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #      Taylor plot.                                                            #
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Plot title.                                                              #
            #------------------------------------------------------------------------------#
            letitre = paste(" Taylor diagram - ",this.desc,"\n",diel.desc[d],sep="")
            cat("       - Taylor","\n")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all formats.                                                  #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.taylor = out[[outform[o]]]$taylor[[diel.key[d]]]
               fichier    = file.path(out.taylor,paste("taylor-",this.vnam,"-",diel.key[d]
                                                      ,".",outform[o],sep=""))
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into 3, and add site and simulation legends at the   #
               # bottom.                                                                   #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,3.0,0))
               layout(mat = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3)),height = c(5.0,1.0))
               #---------------------------------------------------------------------------#




               #----- Legend: the sites. --------------------------------------------------#
               par(mar=c(0.2,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = toupper(sites$iata)
                      , col     = foreground
                      , pt.bg   = foreground
                      , pch     = sites$pch
                      , ncol    = min(4,pretty.box(nsites)$ncol)
                      , title   = expression(bold("Sites"))
                      , pt.cex  = st.cex.med
                      , cex     = 1.1 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#




               #----- Legend: the counts. -------------------------------------------------#
               par(mar=c(0.2,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = c(st.cnt.min,st.cnt.med,st.cnt.max)
                      , col     = foreground
                      , pt.bg   = foreground
                      , pch     = 15
                      , ncol    = 1
                      , title   = expression(bold("Number Obs."))
                      , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                      , cex     = 1.0 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#




               #----- Legend: the simulations. --------------------------------------------#
               par(mar=c(0.2,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = simul$desc
                      , fill    = simul$col
                      , border  = simul$col
                      , ncol    = 2
                      , title   = expression(bold("Structure"))
                      , cex     = 1.0 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Loop over sites.                                                      #
               #---------------------------------------------------------------------------#
               add = FALSE
               for (p in sequence(nsites)){
                  iata = sites$iata[p]

                  #----- Skip the site if there is no data. -------------------------------#
                  ok.iata = (  length(obs.diel[[iata]]) > 0 
                            && any(is.finite(obs.diel[[iata]])) )
                  ok.iata = ok.iata && ( ! is.na(ok.iata))
                  if (ok.iata){
                     #---------------------------------------------------------------------#
                     #     Initialise or update the Taylor plot.                           #
                     #---------------------------------------------------------------------#
                     mytaylor = taylor.plot( obs        = obs.diel[[iata]]
                                           , mod        = mod.diel[[iata]]
                                           , add        = add
                                           , pos.corr   = NA
                                           , pt.col     = simul$col
                                           , pt.bg      = simul$col
                                           , pt.pch     = sites$pch[p]
                                           , pt.cex     = cex.diel [p]
                                           , pt.lwd     = lwd.diel [p]
                                           , obs.col    = foreground
                                           , gamma.col  = sky.mg
                                           , gamma.bg   = background
                                           , sd.col     = grey.fg
                                           , sd.obs.col = yellow.mg
                                           , corr.col   = foreground
                                           , main       = ""
                                           , normalise  = TRUE
                                           )#end taylor.plot
                     add = TRUE
                  }#end if (length(obs.diel[[iata]]) > 0.)
                  #------------------------------------------------------------------------#
               }#end for (p in 1:nsites)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         }#end if (ok.skill.taylor)
         #---------------------------------------------------------------------------------#
      }#end for (d in sequence(ndiel))
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#











      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #      Default simulation, comparing different times of the day.                     #
      #------------------------------------------------------------------------------------#
      cat ("     * Default simulation by time of the day...","\n")
      s         = sim.default
      loop.diel = which(diel.default)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Loop over all sites, normalise data and create the paired vectors.            #
      #------------------------------------------------------------------------------------#
      list.skill    = list()
      percent.skill = array(data=0,dim=c(nsites,ndiel),dimnames=list(sites$iata,diel.key))
      bias.range    = NULL
      sigma.range   = NULL
      for (p in sequence(nsites)){
         #----- Load observation and model. -----------------------------------------------#
         iata  = sites$iata[p]
         obs   = eft[[iata]]
         mod   = res[[iata]]$ans[[simul.key[s]]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Decide which data set to use.                                              #
         #---------------------------------------------------------------------------------#
         ndat    = colSums(is.finite(obs[[this.vnam]]))
         nlyr    = length(ndat)
         slz.use = ifelse(ndat==0,-obs$slz,obs$slz)
         cc      = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
         #---------------------------------------------------------------------------------#


         #----- Initialise list for this site. --------------------------------------------#
         list.skill[[iata]] = list()
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Loop over the diels that will go to the skill/Taylor plots.                  #
         #---------------------------------------------------------------------------------#
         for (d in loop.diel){
            this.diel = diel.key[d]
            #------------------------------------------------------------------------------#
            #      Decide which data set to use.  Either way, use only finite data set,    #
            # and only when the standard deviation is valid and greater than zero          #
            #------------------------------------------------------------------------------#
            if (d %in% diel.fnmean){
               #----- Select fortnightly averages. ----------------------------------------#
               sel                = is.finite(obs[[this.fnmean]][,cc])
               sdev.obs.now       = sd(obs[[this.fnmean]][sel,cc],na.rm=TRUE)
               sel                = sel & sdev.obs.now %>% 0.
               this.obs           = obs[[this.fnmean]][sel,cc]
               this.mod           = mod[[this.fnmean]][sel,cc]
               percent.skill[p,d] = sum(sel)
            }else if(d %in% diel.dmean){
               #----- Select daily averages. ----------------------------------------------#
               sel                = is.finite(obs[[this.dmean]][,cc])
               sdev.obs.now       = sd(obs[[this.dmean]][sel,cc],na.rm=TRUE)
               sel                = sel & sdev.obs.now %>% 0.
               this.obs           = obs[[this.dmean]][sel,cc]
               this.mod           = mod[[this.dmean]][sel,cc]
               percent.skill[p,d] = sum(sel)
            }else{
               #----- Select this diel (or everything for all day). -----------------------#
               d.sel              = obs$diel == d | d %in% diel.all.hrs
               s.sel              = obs$highsun | (! this.sun)
               o.sel              = is.finite(obs[[this.vnam]][,cc])
               nmax               = sum(d.sel & s.sel)
               sel                = d.sel & s.sel & o.sel
               sel                = ifelse(is.na(sel),FALSE,sel)
               sdev.obs.now       = sd(obs[[this.vnam]][sel,cc],na.rm=TRUE)
               sel                = sel & sdev.obs.now %>% 0.
               this.obs           = obs[[this.vnam]][sel,cc]
               this.mod           = mod[[this.vnam]][sel,cc]
               if (sum(d.sel & s.sel) > 0){
                  percent.skill[p,d] = sum(sel)
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Copy the data to the paired list.                                        #
            #------------------------------------------------------------------------------#
            list.skill[[iata]][[this.diel]] = list(obs=this.obs,mod=this.mod)
            #------------------------------------------------------------------------------#


            #----- Find the normalised bias and model standard deviation. -----------------#
            if (percent.skill[p,d] > 0.){
               comp         = res[[iata]]$sim[[simul.key[s]]][[this.vnam]]
               sdev.obs.now = sqrt(comp$obs.moment[cc,d,nseason,2])
               bias.now     = comp$bias [cc,d,nseason] / sdev.obs.now
               sigma.now    = comp$sigma[cc,d,nseason] / sdev.obs.now
               bias.range   = c(bias.range ,bias.now   )
               sigma.range  = c(sigma.range,sigma.now  )
            }#end if (percent.skill[p,d] > 0.)
            #------------------------------------------------------------------------------#
         }#end for (d in loop.diel)
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Plot Taylor and skill plots only if there is anything to plot.                 #
      #------------------------------------------------------------------------------------#
      ok.taylor.skill = (  any(percent.skill > 0)
                        && any(is.finite(bias.range))    && any(is.finite(sigma.range)))
      if (ok.taylor.skill){
         #----- Normalise percentage of data by diel. -------------------------------------#
         diel.scale    = apply(X=percent.skill,MARGIN=2,FUN=max,na.rm=TRUE)
         diel.scale    = array( data     = rep( x    = ifelse( diel.scale == 0
                                                             , 0.
                                                             , 100./diel.scale
                                                             )#end ifelse
                                              , each = nsites)
                              , dim      = dim(percent.skill)
                              , dimnames = dimnames(percent.skill)
                              )#end array
         percent.skill = percent.skill * diel.scale
         #---------------------------------------------------------------------------------#



         #----- Find which combination of place and diel to add to the plots. -------------#
         which.skill = which(x=is.finite(percent.skill) & percent.skill > 0,arr.ind=TRUE)
         nwhich      = nrow(which.skill)
         #---------------------------------------------------------------------------------#




         #---- Fix ranges. ----------------------------------------------------------------#
         xy.range    = 1.04 * max(abs(c(bias.range,sigma.range)),na.rm=TRUE)
         bias.range  = 1.04 * xy.range  * c(-1,1)
         sigma.range = 1.04 * xy.range  * c( 1,0)
         r2.range    = range(1-xy.range^2,1)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Calculate the size for the points in the Skill and Taylor diagrams.         #
         # Make it proportional to the number of points used to evaluate each place.       #
         #---------------------------------------------------------------------------------#
         st.cnt.min = min (percent.skill[percent.skill > 0] , na.rm = TRUE)
         st.cnt.max = max (percent.skill[percent.skill > 0] , na.rm = TRUE)
         st.cnt.med = round(mean(c(st.cnt.min,st.cnt.max)))
         cex.skill  = pmax( st.cex.min, ( st.cex.min + ( st.cex.max    - st.cex.min )
                                                     * ( percent.skill - st.cnt.min )
                                                     / ( st.cnt.max    - st.cnt.min )
                                        )#end cex.skill
                          )#end pmax
         lwd.skill  = pmax( st.lwd.min, ( st.lwd.min + ( st.lwd.max    - st.lwd.min )
                                                     * ( percent.skill - st.cnt.min )
                                                     / ( st.cnt.max    - st.cnt.min )
                                        )#end cex.diel
                         )#end pmax
         cex.skill  = cex.skill + 0. * percent.skill
         lwd.skill  = lwd.skill + 0. * percent.skill
         st.cex.med = mean(c(st.cex.min,st.cex.max))
         #---------------------------------------------------------------------------------#


         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Skill plot.                                                                 #
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Plot title.                                                                 #
         #---------------------------------------------------------------------------------#
         letitre = paste(" Skill diagram - ",this.desc,sep="")
         cat("       - Skill","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.skill = out[[outform[o]]]$skill$default.var
            fichier   = file.path(out.skill,paste("skill-",this.vnam,"-",simul$name[s]
                                                 ,".",outform[o],sep="")
                                 )#end file.path
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
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into 3, and add site and simulation legends at the      #
            # bottom.                                                                      #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,3.0,0))
            layout(mat = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3)),height = c(5.0,1.0))
            #------------------------------------------------------------------------------#




            #----- Legend: the sites. -----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = toupper(sites$iata)
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = sites$pch
                   , ncol    = min(4,pretty.box(nsites)$ncol)
                   , title   = expression(bold("Sites"))
                   , pt.cex  = st.cex.med
                   , cex     = 1.1 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the counts. ----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = sprintf("%.1f",c(st.cnt.min,st.cnt.med,st.cnt.max))
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = 15
                   , ncol    = 1
                   , title   = expression(bold("Rel. Count"))
                   , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the diel. ------------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = diel.desc[diel.default]
                   , fill    = diel.col [diel.default]
                   , border  = diel.col [diel.default]
                   , ncol    = 2
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over valid pairs.                                                   #
            #------------------------------------------------------------------------------#
            myskill     = NULL
            for (u in sequence(nwhich)){
               p         = which.skill[u,1]
               d         = which.skill[u,2]
               iata      = sites$iata[p]
               this.diel = diel.key[d]
               pair.now  = list.skill[[iata]][[this.diel]]


               #---------------------------------------------------------------------------#
               #     Initialise or update the skill plot.                                  #
               #---------------------------------------------------------------------------#
               myskill = skill.plot( obs           = pair.now$obs
                                   , obs.options   = list( col = foreground
                                                         , cex = 2.0
                                                         )#end list
                                   , mod           = pair.now$mod
                                   , mod.options   = list( col = diel.col [d]
                                                         , bg  = diel.col [d]
                                                         , pch = sites$pch[p]
                                                         , cex = cex.skill[p,d]
                                                         , lty = "solid"
                                                         , lwd = lwd.skill[p,d]
                                                         )#end list
                                   , main           = ""
                                   , bias.lim       = bias.range
                                   , r2.lim         = r2.range
                                   , r2.options     = list( col = grid.colour)
                                   , nobias.options = list( col = khaki.mg   )
                                   , rmse.options   = list( col = orange.mg
                                                          , lty = "dotdash"
                                                          , lwd = 1.2
                                                          , bg  = background
                                                          )#end list
                                   , cex.xyzlab     = 1.4
                                   , cex.xyzat      = 1.4
                                   , skill          = myskill
                                   , normalise      = TRUE
                                   , mar            = c(5,4,4,3)+0.1
                                   )#end skill.plot
               #---------------------------------------------------------------------------#
            }#end for (u in sequence(nrow(which.skill)))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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
         }#end for (o in 1:nout) 
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Taylor plot.                                                               #
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Plot title.                                                                 #
         #---------------------------------------------------------------------------------#
         letitre = paste(" Taylor diagram - ",this.desc,sep="")
         cat("       - Taylor","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.taylor = out[[outform[o]]]$taylor$default.var
            fichier    = file.path(out.taylor,paste("taylor-",this.vnam,"-",simul$name[s]
                                                   ,".",outform[o],sep=""))
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
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into 3, and add site and simulation legends at the      #
            # bottom.                                                                      #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,3.0,0))
            layout(mat = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3)),height = c(5.0,1.0))
            #------------------------------------------------------------------------------#




            #----- Legend: the sites. -----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = toupper(sites$iata)
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = sites$pch
                   , ncol    = min(4,pretty.box(nsites)$ncol)
                   , title   = expression(bold("Sites"))
                   , pt.cex  = st.cex.med
                   , cex     = 1.1 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the counts. ----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = sprintf("%.1f",c(st.cnt.min,st.cnt.med,st.cnt.max))
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = 15
                   , ncol    = 1
                   , title   = expression(bold("Rel. Count"))
                   , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the simulations. -----------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = diel.desc[diel.default]
                   , fill    = diel.col [diel.default]
                   , border  = diel.col [diel.default]
                   , ncol    = 2
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over valid pairs.                                                   #
            #------------------------------------------------------------------------------#
            add  = FALSE
            for (u in sequence(nwhich)){
               p         = which.skill[u,1]
               d         = which.skill[u,2]
               iata      = sites$iata[p]
               this.diel = diel.key[d]
               pair.now  = list.skill[[iata]][[this.diel]]

               #---------------------------------------------------------------------------#
               #     Initialise or update the Taylor plot.                                 #
               #---------------------------------------------------------------------------#
               mytaylor = taylor.plot( obs        = pair.now$obs
                                     , mod        = pair.now$mod
                                     , add        = add
                                     , pos.corr   = NA
                                     , pt.col     = diel.col [d]
                                     , pt.bg      = diel.col [d]
                                     , pt.pch     = sites$pch[p]
                                     , pt.cex     = cex.skill[p,d]
                                     , pt.lwd     = lwd.skill[p,d]
                                     , obs.col    = foreground
                                     , gamma.col  = sky.mg
                                     , gamma.bg   = background
                                     , sd.col     = grey.fg
                                     , sd.obs.col = yellow.mg
                                     , corr.col   = foreground
                                     , main       = ""
                                     , normalise  = TRUE
                                     )#end taylor.plot
               add = TRUE
               #---------------------------------------------------------------------------#
            }#end for (u in sequence(nwhich))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end if (ok.taylor.skill)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #       Default simulation, all variables in the same plot, colours are sites, point    #
   # formats are variables.                                                                #
   #---------------------------------------------------------------------------------------#
   cat ("     * Default simulation by time of the day...","\n")
   s         = sim.default
   loop.diel = which(diel.default)
   v.inc     = c("gpp","nep","reco","rshortup","rlongup","ustar","hflxca","wflxca")
   v.use     = match(v.inc,compvar.key)
   nvuse     = length(v.use)

   for (d in sequence(ndiel)){
      this.diel = diel.key[d]

      #------------------------------------------------------------------------------------#
      #      Loop over all sites, normalise data and create the paired vectors.            #
      #------------------------------------------------------------------------------------#
      list.skill    = list()
      percent.skill = array( data      = 0
                           , dim       = c(ncompvar,nsites)
                           , dimnames  = list(compvar.key,sites$iata)
                           )#end array
      bias.range    = NULL
      sigma.range   = NULL
      for (v in sequence(ncompvar)){
         #----- Copy the variable information. --------------------------------------------#
         this.vnam     = compvar[[v]]$vnam
         this.dmean    = paste("dmean"   ,this.vnam,sep=".")
         this.fnmean   = paste("fnmean"  ,this.vnam,sep=".")
         this.measured = paste("measured",this.vnam,sep=".")
         this.desc     = compvar[[v]]$desc
         this.unit     = compvar[[v]]$unit
         this.sun      = compvar[[v]]$sunvar
         this.soilvar  = compvar[[v]]$soilvar
         #---------------------------------------------------------------------------------#


         #----- Initialise list for this variable. ----------------------------------------#
         list.skill[[this.vnam]] = list()
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Site loop.                                                                 #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            #----- Load observation and model. --------------------------------------------#
            iata  = sites$iata[p]
            obs   = eft[[iata]]
            mod   = res[[iata]]$ans[[simul.key[s]]]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Decide which data set to use.                                           #
            #------------------------------------------------------------------------------#
            ndat    = colSums(is.finite(obs[[this.vnam]]))
            nlyr    = length(ndat)
            slz.use = ifelse(ndat==0,-obs$slz,obs$slz)
            cc      = pmin(nlyr,which.min(abs(slz.use-slz.reference)))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Decide which data set to use.  Either way, use only finite data set,    #
            # and only when the standard deviation is valid and greater than zero          #
            #------------------------------------------------------------------------------#
            if (d %in% diel.fnmean){
               #----- Select fortnightly averages. ----------------------------------------#
               sel                = ( is.finite(obs[[this.fnmean]][,cc])
                                    & this.vnam %in% v.inc )
               sdev.obs.now       = sd(obs[[this.fnmean]][sel,cc],na.rm=TRUE)
               sel                = sel & sdev.obs.now %>% 0.
               this.obs           = obs[[this.fnmean]][sel,cc]
               this.mod           = mod[[this.fnmean]][sel,cc]
               percent.skill[v,p] = sum(sel)
            }else if(d %in% diel.dmean){
               #----- Select daily averages. ----------------------------------------------#
               sel                = ( is.finite(obs[[this.dmean]][,cc])
                                    & this.vnam %in% v.inc )
               sdev.obs.now       = sd(obs[[this.dmean]][sel,cc],na.rm=TRUE)
               sel                = sel & sdev.obs.now %>% 0.
               this.obs           = obs[[this.dmean]][sel,cc]
               this.mod           = mod[[this.dmean]][sel,cc]
               percent.skill[v,p] = sum(sel)
            }else{
               #----- Select this diel (or everything for all day). -----------------------#
               d.sel              = obs$diel == d | d %in% diel.all.hrs
               s.sel              = obs$highsun | (! this.sun)
               o.sel              = is.finite(obs[[this.vnam]][,cc]) & this.vnam %in% v.inc
               nmax               = sum(d.sel & s.sel)
               sel                = d.sel & s.sel & o.sel
               sel                = ifelse(is.na(sel),FALSE,sel)
               sdev.obs.now       = sd(obs[[this.vnam]][sel,cc],na.rm=TRUE)
               sel                = sel & sdev.obs.now %>% 0.
               this.obs           = obs[[this.vnam]][sel,cc]
               this.mod           = mod[[this.vnam]][sel,cc]
               if (sum(d.sel & s.sel) > 0){
                  percent.skill[v,p] = sum(sel)
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Copy the data to the paired list.                                        #
            #------------------------------------------------------------------------------#
            list.skill[[this.vnam]][[iata]] = list(obs=this.obs,mod=this.mod)
            #------------------------------------------------------------------------------#



            #----- Find the normalised bias and model standard deviation. -----------------#
            if (percent.skill[v,p] > 0.){
               comp         = res[[iata]]$sim[[simul.key[s]]][[this.vnam]]
               sdev.obs.now = sqrt(comp$obs.moment[cc,d,nseason,2])
               bias.now     = comp$bias [cc,d,nseason] / sdev.obs.now
               sigma.now    = comp$sigma[cc,d,nseason] / sdev.obs.now
               bias.range   = c(bias.range ,bias.now   )
               sigma.range  = c(sigma.range,sigma.now  )
            }#end if (percent.skill[v,p] > 0.)
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#
      }#end for (v in v.use)
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Plot Taylor and skill plots only if there is anything to plot.                 #
      #------------------------------------------------------------------------------------#
      ok.taylor.skill = (  any(percent.skill > 0)
                        && any(is.finite(bias.range))    && any(is.finite(sigma.range)))
      if (ok.taylor.skill){
         #----- Normalise percentage of data by diel. -------------------------------------#
         diel.scale    = apply(X=percent.skill,MARGIN=1,FUN=max,na.rm=TRUE)
         diel.scale    = array( data     = rep( x     = ifelse( diel.scale == 0
                                                              , 0.
                                                              , 100./diel.scale
                                                              )#end ifelse
                                              , times = nsites)
                              , dim      = dim(percent.skill)
                              , dimnames = dimnames(percent.skill)
                              )#end array
         percent.skill = percent.skill * diel.scale
         #---------------------------------------------------------------------------------#



         #----- Find which combination of place and diel to add to the plots. -------------#
         which.skill = which(x=is.finite(percent.skill) & percent.skill > 0,arr.ind=TRUE)
         nwhich      = nrow(which.skill)
         #---------------------------------------------------------------------------------#




         #---- Fix ranges. ----------------------------------------------------------------#
         xy.range    = 1.04 * max(abs(c(bias.range,sigma.range)),na.rm=TRUE)
         bias.range  = 1.04 * xy.range  * c(-1,1)
         sigma.range = 1.04 * xy.range  * c( 1,0)
         r2.range    = range(1-xy.range^2,1)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Calculate the size for the points in the Skill and Taylor diagrams.         #
         # Make it proportional to the number of points used to evaluate each place.       #
         #---------------------------------------------------------------------------------#
         st.cnt.min = min (percent.skill[percent.skill > 0] , na.rm = TRUE)
         st.cnt.max = max (percent.skill[percent.skill > 0] , na.rm = TRUE)
         st.cnt.med = round(mean(c(st.cnt.min,st.cnt.max)))
         cex.skill  = pmax( st.cex.min, ( st.cex.min + ( st.cex.max    - st.cex.min )
                                                     * ( percent.skill - st.cnt.min )
                                                     / ( st.cnt.max    - st.cnt.min )
                                        )#end cex.skill
                          )#end pmax
         lwd.skill  = pmax( st.lwd.min, ( st.lwd.min + ( st.lwd.max    - st.lwd.min )
                                                     * ( percent.skill - st.cnt.min )
                                                     / ( st.cnt.max    - st.cnt.min )
                                        )#end cex.diel
                         )#end pmax
         cex.skill  = cex.skill + 0. * percent.skill
         lwd.skill  = lwd.skill + 0. * percent.skill
         st.cex.med = mean(c(st.cex.min,st.cex.max))
         #---------------------------------------------------------------------------------#


         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Skill plot.                                                                 #
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Plot title.                                                                 #
         #---------------------------------------------------------------------------------#
         letitre = paste(" Skill diagram - ",diel.desc[d],sep="")
         cat("       - Skill","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.skill = out[[outform[o]]]$skill$default.site
            fichier   = file.path(out.skill,paste("skill-",diel.key[d],"-",simul$name[s]
                                                 ,".",outform[o],sep="")
                                 )#end file.path
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
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into 3, and add site and simulation legends at the      #
            # bottom.                                                                      #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,3.0,0))
            layout(mat = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3)),height = c(5.0,1.0))
            #------------------------------------------------------------------------------#




            #----- Legend: the variables. -------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = compvar.sym[v.use]
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = compvar.pch[v.use]
                   , ncol    = min(4,pretty.box(nvuse)$ncol)
                   , title   = expression(bold("Sites"))
                   , pt.cex  = st.cex.med
                   , cex     = 1.1 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the counts. ----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = sprintf("%.1f",c(st.cnt.min,st.cnt.med,st.cnt.max))
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = 15
                   , ncol    = 1
                   , title   = expression(bold("Rel. Count"))
                   , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the sites. -----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = paste(sites$desc," (",toupper(sites$iata),")",sep="")
                   , fill    = sites$col
                   , border  = sites$col
                   , ncol    = 2
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over valid pairs.                                                   #
            #------------------------------------------------------------------------------#
            myskill     = NULL
            for (u in sequence(nwhich)){
               v         = which.skill[u,1]
               p         = which.skill[u,2]
               this.vnam = compvar[[v]]$vnam
               iata      = sites$iata[p]
               pair.now  = list.skill[[this.vnam]][[iata]]


               #---------------------------------------------------------------------------#
               #     Initialise or update the skill plot.                                  #
               #---------------------------------------------------------------------------#
               myskill = skill.plot( obs           = pair.now$obs
                                   , obs.options   = list( col = foreground
                                                         , cex = 2.0
                                                         )#end list
                                   , mod           = pair.now$mod
                                   , mod.options   = list( col = sites$col[p]
                                                         , bg  = sites$col[p]
                                                         , pch = compvar[[v]]$pch
                                                         , cex = cex.skill[v,p]
                                                         , lty = "solid"
                                                         , lwd = lwd.skill[v,p]
                                                         )#end list
                                   , main           = ""
                                   , bias.lim       = bias.range
                                   , r2.lim         = r2.range
                                   , r2.options     = list( col = grid.colour)
                                   , nobias.options = list( col = khaki.mg   )
                                   , rmse.options   = list( col = orange.mg
                                                          , lty = "dotdash"
                                                          , lwd = 1.2
                                                          , bg  = background
                                                          )#end list
                                   , cex.xyzlab     = 1.4
                                   , cex.xyzat      = 1.4
                                   , skill          = myskill
                                   , normalise      = TRUE
                                   , mar            = c(5,4,4,3)+0.1
                                   )#end skill.plot
               #---------------------------------------------------------------------------#
            }#end for (u in sequence(nrow(which.skill)))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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
         }#end for (o in 1:nout) 
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Taylor plot.                                                               #
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Plot title.                                                                 #
         #---------------------------------------------------------------------------------#
         letitre = paste(" Taylor diagram - ",diel.desc[d],sep="")
         cat("       - Taylor","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.taylor = out[[outform[o]]]$taylor$default.site
            fichier    = file.path(out.taylor,paste("taylor-",diel.key[d],"-",simul$name[s]
                                                   ,".",outform[o],sep=""))
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
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into 3, and add site and simulation legends at the      #
            # bottom.                                                                      #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,3.0,0))
            layout(mat = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3)),height = c(5.0,1.0))
            #------------------------------------------------------------------------------#




            #----- Legend: the variables. -------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = compvar.sym[v.use]
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = compvar.pch[v.use]
                   , ncol    = min(4,pretty.box(nvuse)$ncol)
                   , title   = expression(bold("Sites"))
                   , pt.cex  = st.cex.med
                   , cex     = 1.1 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the counts. ----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = sprintf("%.1f",c(st.cnt.min,st.cnt.med,st.cnt.max))
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = 15
                   , ncol    = 1
                   , title   = expression(bold("Rel. Count"))
                   , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Legend: the sites. -----------------------------------------------------#
            par(mar=c(0.2,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = paste(sites$desc," (",toupper(sites$iata),")",sep="")
                   , fill    = sites$col
                   , border  = sites$col
                   , ncol    = 2
                   , cex     = 1.0 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over valid pairs.                                                   #
            #------------------------------------------------------------------------------#
            add  = FALSE
            for (u in sequence(nwhich)){
               v         = which.skill[u,1]
               p         = which.skill[u,2]
               this.vnam = compvar[[v]]$vnam
               iata      = sites$iata[p]
               pair.now  = list.skill[[this.vnam]][[iata]]

               #---------------------------------------------------------------------------#
               #     Initialise or update the Taylor plot.                                 #
               #---------------------------------------------------------------------------#
               mytaylor = taylor.plot( obs        = pair.now$obs
                                     , mod        = pair.now$mod
                                     , add        = add
                                     , pos.corr   = NA
                                     , pt.col     = sites$col[p]
                                     , pt.bg      = sites$col[p]
                                     , pt.pch     = compvar[[v]]$pch
                                     , pt.cex     = cex.skill[v,p]
                                     , pt.lwd     = lwd.skill[v,p]
                                     , obs.col    = foreground
                                     , gamma.col  = sky.mg
                                     , gamma.bg   = background
                                     , sd.col     = grey.fg
                                     , sd.obs.col = yellow.mg
                                     , corr.col   = foreground
                                     , main       = ""
                                     , normalise  = TRUE
                                     )#end taylor.plot
               add = TRUE
               #---------------------------------------------------------------------------#
            }#end for (u in sequence(nwhich))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
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
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end if (ok.taylor.skill)
      #------------------------------------------------------------------------------------#
   }#end for (d in loop.diel)
   #---------------------------------------------------------------------------------------#
}#end if (plot.skill.taylor)
#------------------------------------------------------------------------------------------#












#------------------------------------------------------------------------------------------#
#      Plot the skill as a function of depth.                                              #
#------------------------------------------------------------------------------------------#
if (plot.soil.skill){
   cat (" + Plot skill plots for soil variables...","\n")
   for (v in sequence(ncompvar)){
      #----- Copy the variable information. -----------------------------------------------#
      this.vnam     = compvar[[v]]$vnam
      this.dmean    = paste("dmean"   ,this.vnam,sep=".")
      this.fnmean   = paste("fnmean"  ,this.vnam,sep=".")
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = compvar[[v]]$desc
      this.unit     = compvar[[v]]$unit
      this.sun      = compvar[[v]]$sunvar
      this.soilvar  = compvar[[v]]$soilvar
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Skip if this isn't a soil variable.                                            #
      #------------------------------------------------------------------------------------#
      if (this.soilvar){
         cat("   - ",this.desc,"...","\n")




         #---------------------------------------------------------------------------------#
         #      Loop over all parts of the day.                                            #
         #---------------------------------------------------------------------------------#
         for (d in sequence(ndiel)){
            cat("     * ",diel.desc[d],"...","\n")


            #------------------------------------------------------------------------------#
            #      Loop over all sites, and get every layer that has data.                 #
            #------------------------------------------------------------------------------#
            obs.diel    = list()
            mod.diel    = list()
            cnt.diel    = list()
            bias.iata   = mapply(FUN=numeric,rep(0,times=nsites),SIMPLIFY=FALSE)
            sigma.iata  = mapply(FUN=numeric,rep(0,times=nsites),SIMPLIFY=FALSE)
            bias.simul  = mapply(FUN=numeric,rep(0,times=nsimul),SIMPLIFY=FALSE)
            sigma.simul = mapply(FUN=numeric,rep(0,times=nsimul),SIMPLIFY=FALSE)
            names(bias.iata)  = names(sigma.iata)  = sites$iata
            names(bias.simul) = names(sigma.simul) = simul$name
            for (p in sequence(nsites)){
               #----- Grab data. ----------------------------------------------------------#
               iata              = sites$iata[p]
               obs               = eft[[iata]]
               nwhen             = obs$nwhen
               obs.diel[[iata]]  = list()
               mod.diel[[iata]]  = list()
               cnt.diel[[iata]]  = list()
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Decide which data set to use.                                        #
               #---------------------------------------------------------------------------#
               ndat     = colSums(is.finite(obs[[this.vnam]]))
               nlyr     = length(ndat)
               loop.nzg = sequence(nlyr)
               for (cc in loop.nzg){
                  if (d %in% diel.fnmean){
                     this.obs = obs[[this.fnmean]][,cc]
                     sel      = is.finite(this.obs)
                     this.obs = this.obs[sel]
                  }else if (d %in% diel.dmean){
                     this.obs = obs[[this.dmean]][,cc]
                     sel      = is.finite(this.obs)
                     this.obs = this.obs[sel]
                  }else{
                     #----- Select this diel (or everything for all day). -----------------#
                     d.sel    = (obs$diel == d | d %in% diel.all.hrs)
                     s.sel    = obs$highsun | (! this.sun)
                     o.sel    = is.finite(obs[[this.vnam]][,cc])
                     sel      = d.sel & s.sel & o.sel
                     sel      = ifelse(is.na(sel),FALSE,sel)
                     this.obs = obs[[this.vnam]][sel,cc]
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Find the standard deviation of this observation.  Skip the site if #
                  # everything is zero.                                                    #
                  #------------------------------------------------------------------------#
                  comp         = res[[iata]]$sim[[simul.key[1]]][[this.vnam]]
                  sdev.obs.now = sqrt(comp$obs.moment[cc,d,nseason,2])
                  sel          = sel & sdev.obs.now %>% 0
                  n.sel        = sum(sel)
                  #------------------------------------------------------------------------#



                  #----- Copy the observed data. ------------------------------------------#
                  obs.diel[[iata]][[cc]] = this.obs
                  mod.diel[[iata]][[cc]] = matrix( ncol     = nsimul
                                                 , nrow     = n.sel
                                                 , dimnames = list(NULL,simul.key)
                                                 )#end matrix
                  cnt.diel[[iata]][[cc]] = sum(is.finite(this.obs))
                  #------------------------------------------------------------------------#



                  #----- Copy the modelled data, and update ranges. -----------------------#
                  if (any(sel)){
                     for (s in sequence(nsimul)){


                        #------------------------------------------------------------------#
                        #      Decide which variable to use depending on the diel.         #
                        #------------------------------------------------------------------#
                        mod  = res[[iata]]$ans[[simul.key[s]]]
                        if (d %in% diel.fnmean){
                           this.mod = mod[[this.fnmean]][sel,cc]
                        }else if (d %in% diel.dmean){
                           this.mod = mod[[this.dmean]][sel,cc]
                        }else{
                           this.mod = mod[[this.vnam]][sel,cc]
                        }#end if
                        mod.diel[[iata]][[cc]][,s] = this.mod
                        #------------------------------------------------------------------#



                        #----- Find the normalised bias and model standard deviation. -----#
                        comp                = res[[iata]]$sim[[simul.key[s]]][[this.vnam]]
                        sdev.obs.now        = sqrt(comp$obs.moment[cc,d,nseason,2])
                        bias.now            = comp$bias [cc,d,nseason] / sdev.obs.now
                        sigma.now           = comp$sigma[cc,d,nseason] / sdev.obs.now
                        bias.iata  [[iata]] = c(bias.iata  [[iata]],bias.now )
                        sigma.iata [[iata]] = c(sigma.iata [[iata]],sigma.now)
                        bias.simul [[s   ]] = c(bias.simul [[s   ]],bias.now )
                        sigma.simul[[s   ]] = c(sigma.simul[[s   ]],sigma.now)
                        #------------------------------------------------------------------#
                     }#end for (s in sequence(nsimul))
                     #---------------------------------------------------------------------#
                  }#end if (any(sel))
                  #------------------------------------------------------------------------#
               }#end for (cc in loop.nzg)
               #---------------------------------------------------------------------------#









               #---------------------------------------------------------------------------#
               #     Plot Taylor and skill plots only if there is anything to plot.        #
               #---------------------------------------------------------------------------#
               bias.range  = bias.iata [[iata]]
               sigma.range = sigma.iata[[iata]]
               ok.taylor.skill = (  length(unlist(obs.diel[[iata]])) > 0
                                 && any(cnt.diel[[iata]] > 0)
                                 && any(is.finite(bias.range ))
                                 && any(is.finite(sigma.range)))
               if (ok.taylor.skill){
                  #---- Fix ranges. -------------------------------------------------------#
                  xy.range    = 1.04 * max(abs(c(bias.range,sigma.range)),na.rm=TRUE)
                  bias.range  = 1.04 * xy.range  * c(-1,1)
                  sigma.range = 1.04 * xy.range  * c( 1,0)
                  r2.range    = range(1-xy.range^2,1)
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Calculate the size for the points in the Skill and Taylor          #
                  # diagrams.  Make it proportional to the number of points used to        #
                  # evaluate each place.                                                   #
                  #------------------------------------------------------------------------#
                  cnt.now    = cnt.diel[[iata]]
                  st.cnt.min = min(unlist(cnt.now),na.rm=TRUE)
                  st.cnt.max = max(unlist(cnt.now),na.rm=TRUE)
                  st.cnt.med = round(mean(c(st.cnt.min,st.cnt.max)))
                  st.cex.med = mean(c(st.cex.min,st.cex.max))
                  cex.diel   = list()
                  lwd.diel   = list()
                  for (cc in seq_along(cnt.now)){
                     cex.diel[[cc]] = pmax( st.cex.min
                                          , ( st.cex.min 
                                            + ( st.cex.max     - st.cex.min )
                                            * ( cnt.now[[cc]]  - st.cnt.min )
                                            / ( st.cnt.max     - st.cnt.min )
                                            )#end cex.diel
                                          )#end pmax
                     lwd.diel[[cc]] = pmax( st.lwd.min
                                          , ( st.lwd.min
                                            + ( st.lwd.max     - st.lwd.min )
                                            * ( cnt.now[[cc]]  - st.cnt.min )
                                            / ( st.cnt.max     - st.cnt.min )
                                            )#end cex.diel
                                          )#end pmax
                  }#end for
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot title.                                                        #
                  #------------------------------------------------------------------------#
                  letitre = paste(" Skill diagram - ",this.desc,"\n"
                                 ,sites$desc[p]," - ",diel.desc[d],sep="")
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Find the legend for layers.                                       #
                  #------------------------------------------------------------------------#
                  slz.desc = paste(sprintf("%.0f",-100*obs$slz),"cm")
                  slz.col  = get(slz.cscheme)(n=nlyr)
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Loop over all formats.                                            #
                  #------------------------------------------------------------------------#
                  for (o in sequence(nout)){
                     #----- Make the file name. -------------------------------------------#
                     out.skill = out[[outform[o]]]$soil.skill$variables
                     fichier   = file.path(out.skill
                                          ,paste("skill-",this.vnam,"-",diel.key[d],"-"
                                          ,sites$iata[p],".",outform[o],sep="")
                                          )#end file.path
                     if (outform[o] == "x11"){
                        X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=wsize$width*depth
                           ,height=wsize$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=wsize$width,height=wsize$height
                                  ,pointsize=ptsz,paper=size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE,width=wsize$width
                           ,height=wsize$height,pointsize=ptsz,paper=size$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Make the legends at the bottom.                                 #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     par.orig = par(no.readonly = TRUE)
                     mar.orig = par.orig$mar
                     par(oma = c(0.2,3,3.0,0))
                     layout( mat     = rbind(lo.simul$mat.off2,c(1,2))
                           , heights = c(rep(5.0/lo.simul$nrow,lo.simul$nrow),1.0)
                           )#end layout
                     #---------------------------------------------------------------------#




                     #----- Legend: the counts. -------------------------------------------#
                     par(mar=c(0.2,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = c(st.cnt.min,st.cnt.med,st.cnt.max)
                            , col     = foreground
                            , pt.bg   = foreground
                            , pch     = 15
                            , pt.lwd  = 3
                            , ncol    = 3
                            , title   = expression(bold("Number Obs."))
                            , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                            , cex     = 1.0 * cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#




                     #----- Legend: the soil layers. --------------------------------------#
                     par(mar=c(0.2,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = slz.desc
                            , col     = slz.col
                            , pt.lwd  = 3
                            , pch     = 4
                            , ncol    = min(3,pretty.box(length(slz.desc))$ncol)
                            , title   = expression(bold("Soil layer"))
                            , cex     = 0.8 * cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #     Loop over simulations and layers.                               #
                     #---------------------------------------------------------------------#
                     for (s in sequence(nsimul)){
                        myskill = NULL
                        for (cc in seq_along(obs.diel[[iata]])){


                           #---------------------------------------------------------------#
                           #     Make sure that something is going to be on the objects.   #
                           #---------------------------------------------------------------#
                           if (cnt.diel[[iata]][[cc]] != 0){
                              this.obs = obs.diel[[iata]][[cc]]
                              this.mod = mod.diel[[iata]][[cc]][,s]
                           }else{
                              this.obs = rep(NA,times=4)
                              this.mod = rep(NA,times=4)
                           }#end if
                           #---------------------------------------------------------------#


                           #---------------------------------------------------------------#
                           #     Initialise or update the skill plot.                      #
                           #---------------------------------------------------------------#
                           myskill = skill.plot( obs           = this.obs
                                               , obs.options   = list( col = foreground
                                                                     , cex = 2.0
                                                                     )#end list
                                               , mod           = this.mod
                                               , mod.options   = list( col = slz.col[cc]
                                                                     , bg  = slz.col[cc]
                                                                     , pch = 4
                                                                     , cex = cex.diel[[cc]]
                                                                     , lty = "solid"
                                                                     , lwd = lwd.diel[[cc]]
                                                                     )#end list
                                               , main           = simul$desc[s]
                                               , bias.lim       = bias.range
                                               , r2.lim         = r2.range
                                               , r2.options     = list( col = grid.colour)
                                               , nobias.options = list( col = khaki.mg   )
                                               , rmse.options   = list( col = orange.mg
                                                                      , lty = "dotdash"
                                                                      , lwd = 1.2
                                                                      , bg  = background
                                                                      )#end list
                                               , cex.xyzlab     = 1.0
                                               , cex.xyzat      = 1.0
                                               , skill          = myskill
                                               , normalise      = TRUE
                                               , mar            = c(3,3,3,3)+0.1
                                               )#end skill.plot
                           #---------------------------------------------------------------#
                        }#end for (cc in seq_along(obs.diel[[iata]]))
                        #------------------------------------------------------------------#
                     }#end for (s in sequence(nsimul))
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the global title.                                          #
                     #---------------------------------------------------------------------#
                     par(las=0)
                     mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     dummy = clean.tmp()
                     #---------------------------------------------------------------------#
                  }#end for (o in 1:nout)
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
               }#end if (ok.taylor.skill)
               #---------------------------------------------------------------------------#
            }#end for (p in sequence(nsites))
            #------------------------------------------------------------------------------#








            #------------------------------------------------------------------------------#
            #     Plot default simulation.                                                 #
            #------------------------------------------------------------------------------#
            s               = sim.default
            bias.range      = bias.simul [[s]]
            sigma.range     = sigma.simul[[s]]
            #----- Check whether there is anything to plot, and grab all soil layers. -----#
            ok.taylor.skill = FALSE
            slz.all         = NULL
            for (p in sequence(nsites)){
               iata            = sites$iata[p]
               obs             = eft[[iata]]
               ok.taylor.skill = (   ok.taylor.skill
                                 || (  length(unlist(obs.diel[[iata]])) > 0
                                    && any(cnt.diel[[iata]] > 0)
                                    && any(is.finite(bias.range ))
                                    && any(is.finite(sigma.range)))
                                 )
               slz.all = c(slz.all,pmin(-0.05,obs$slz))
            }#end for (p in sequence(nsites))
            slz.all = sort(unique(slz.all))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot Taylor and skill plots only if there is anything to plot.           #
            #------------------------------------------------------------------------------#
            if (ok.taylor.skill){
               #---- Fix ranges. ----------------------------------------------------------#
               xy.range    = 1.04 * max(abs(c(bias.range,sigma.range)),na.rm=TRUE)
               bias.range  = 1.04 * xy.range  * c(-1,1)
               sigma.range = 1.04 * xy.range  * c( 1,0)
               r2.range    = range(1-xy.range^2,1)
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Calculate the size for the points in the Skill and Taylor             #
               # diagrams.  Make it proportional to the number of points used to           #
               # evaluate each place.                                                      #
               #---------------------------------------------------------------------------#
               st.cnt.min = +Inf
               st.cnt.max = -Inf
               for (p in sequence(nsites)){
                  iata       = sites$iata[p]
                  cnt.now    = cnt.diel[[iata]]
                  cnt.ref    = unlist(cnt.now)
                  cnt.ref    = ifelse(cnt.ref==0,NA,cnt.ref)
                  st.cnt.min = min(c(st.cnt.min,cnt.ref),na.rm=TRUE)
                  st.cnt.max = max(c(st.cnt.max,cnt.ref),na.rm=TRUE)
               }#end for
               st.cnt.med = round(mean(c(st.cnt.min,st.cnt.max)))
               st.cex.med = mean(c(st.cex.min,st.cex.max))
               cex.diel   = list()
               lwd.diel   = list()
               for (p in sequence(nsites)){
                  iata = sites$iata[p]
                  cex.diel[[iata]]   = list()
                  lwd.diel[[iata]]   = list()
                  cnt.now    = cnt.diel[[iata]]
                  for (cc in seq_along(cnt.now)){
                     cex.diel[[iata]][[cc]] = pmax( st.cex.min
                                                  , ( st.cex.min 
                                                    + ( st.cex.max     - st.cex.min )
                                                    * ( cnt.now[[cc]]  - st.cnt.min )
                                                    / ( st.cnt.max     - st.cnt.min )
                                                    )#end cex.diel
                                                  )#end pmax
                     lwd.diel[[iata]][[cc]] = pmax( st.lwd.min
                                                  , ( st.lwd.min
                                                    + ( st.lwd.max     - st.lwd.min )
                                                    * ( cnt.now[[cc]]  - st.cnt.min )
                                                    / ( st.cnt.max     - st.cnt.min )
                                                    )#end cex.diel
                                                  )#end pmax
                  }#end for (cc in seq_along(cnt.now))
                  #------------------------------------------------------------------------#
               }#end for (p in sequence(nsites))
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Plot title.                                                           #
               #---------------------------------------------------------------------------#
               letitre = paste(" Skill diagram - ",this.desc,"\n",diel.desc[d],sep="")
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Find the legend for layers.                                          #
               #---------------------------------------------------------------------------#
               slz.desc = paste(sprintf("%.0f",-100*slz.all),"cm")
               slz.col  = get(slz.cscheme)(n=length(slz.all))
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Loop over all formats.                                               #
               #---------------------------------------------------------------------------#
               for (o in sequence(nout)){
                  #----- Make the file name. ----------------------------------------------#
                  out.skill = out[[outform[o]]]$soil.skill$default
                  fichier   = file.path(out.skill
                                       ,paste("skill-",this.vnam,"-",diel.key[d]
                                       ,".",outform[o],sep="")
                                       )#end file.path
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
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Make three legends at the bottom.                                  #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  mar.orig = par.orig$mar
                  par(oma = c(0.2,3,3.0,0))
                  layout(mat     = rbind(c(4,4,4,4,4,4,4),c(1,1,2,3,3,3,3))
                       , heights = c(5.0,1.0)
                       )#end layout
                  #------------------------------------------------------------------------#




                  #----- Legend: the sites. -----------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = toupper(sites$iata)
                         , col     = foreground
                         , pt.bg   = foreground
                         , pch     = sites$pch
                         , ncol    = min(3,pretty.box(nsites)$ncol)
                         , title   = expression(bold("Sites"))
                         , pt.cex  = st.cex.med
                         , cex     = 1.1 * cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#




                  #----- Legend: the counts. ----------------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = c(st.cnt.min,st.cnt.med,st.cnt.max)
                         , col     = foreground
                         , pt.bg   = foreground
                         , pch     = 15
                         , pt.lwd  = 3
                         , ncol    = 1
                         , title   = expression(bold("Number Obs."))
                         , pt.cex  = c(st.cex.min,st.cex.med,st.cex.max)
                         , cex     = 1.0 * cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#




                  #----- Legend: the soil layers. -----------------------------------------#
                  par(mar=c(0.2,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  legend ( x       = "bottom"
                         , inset   = 0.0
                         , legend  = slz.desc
                         , fill    = slz.col
                         , border  = slz.col
                         , ncol    = 5
                         , title   = expression(bold("Soil layer"))
                         , cex     = 0.8 * cex.ptsz
                         , xpd     = TRUE
                         )#end legend
                  #------------------------------------------------------------------------#





                  #------------------------------------------------------------------------#
                  #     Loop over simulations and layers.                                  #
                  #------------------------------------------------------------------------#
                  myskill = NULL
                  for (p in sequence(nsites)){
                     iata    = sites$iata[p]
                     slz.idx = match(pmin(-0.05,eft[[iata]]$slz),slz.all)


                     for (cc in seq_along(obs.diel[[iata]])){
                        #------------------------------------------------------------------#
                        #     Make sure that something is going to be on the objects.      #
                        #------------------------------------------------------------------#
                        if (cnt.diel[[iata]][[cc]] != 0){
                           this.obs = obs.diel[[iata]][[cc]]
                           this.mod = mod.diel[[iata]][[cc]][,s]
                        }else{
                           this.obs = rep(NA,times=4)
                           this.mod = rep(NA,times=4)
                        }#end if
                        this.cex = cex.diel[[iata]][[cc]]
                        this.lwd = lwd.diel[[iata]][[cc]]
                        this.col = slz.col[slz.idx[cc]]
                        #------------------------------------------------------------------#


                        #------------------------------------------------------------------#
                        #     Initialise or update the skill plot.                         #
                        #------------------------------------------------------------------#
                        myskill = skill.plot( obs           = this.obs
                                            , obs.options   = list( col = foreground
                                                                  , cex = 2.0
                                                                  )#end list
                                            , mod           = this.mod
                                            , mod.options   = list( col = this.col
                                                                  , bg  = this.col
                                                                  , pch = sites$pch[p]
                                                                  , cex = this.cex
                                                                  , lty = "solid"
                                                                  , lwd = this.lwd
                                                                  )#end list
                                            , main           = ""
                                            , bias.lim       = bias.range
                                            , r2.lim         = r2.range
                                            , r2.options     = list( col = grey.fg    )
                                            , nobias.options = list( col = grey.fg    )
                                            , rmse.options   = list( col = grey.mg
                                                                   , lty = "dotdash"
                                                                   , lwd = 1.2
                                                                   , bg  = background
                                                                   , cex = 0.8
                                                                   )#end list
                                            , cex.xyzlab     = 1.3
                                            , cex.xyzat      = 1.3
                                            , skill          = myskill
                                            , normalise      = TRUE
                                            , mar            = c(5,4,4,3)+0.1
                                            )#end skill.plot
                        #------------------------------------------------------------------#
                     }#end for (cc in seq_along(obs.diel[[iata]]))
                     #---------------------------------------------------------------------#
                  }#end for (p in sequence(nsites))
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the global title.                                             #
                  #------------------------------------------------------------------------#
                  par(las=0)
                  mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  dummy = clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            }#end if (ok.taylor.skill)
            #------------------------------------------------------------------------------#


         }#end for (d in loop.diel)
         #---------------------------------------------------------------------------------#
      }#end if (this.soilvar)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar)
   #---------------------------------------------------------------------------------------#
}#end if (plot.soil.skill)
#------------------------------------------------------------------------------------------#












#------------------------------------------------------------------------------------------#
#     Make a summary table with all variables and results.                                 #
#------------------------------------------------------------------------------------------#
summ.table = list()
if (make.summ.table){
   cat (" + Make summary table with all observations and sites...","\n")

   #----- Initialise summary table. -------------------------------------------------------#
   template.table    = data.frame( n                = rep(0.           ,times=ncompvar)
                                 , df               = rep(0.           ,times=ncompvar)
                                 , bias             = rep(0.           ,times=ncompvar)
                                 , sigma            = rep(0.           ,times=ncompvar)
                                 , ks.stat          = rep(0,           ,times=ncompvar)
                                 , best.rmse        = rep(Inf          ,times=ncompvar)
                                 , best.iata        = rep(NA_character_,times=ncompvar)
                                 , best.n           = rep(NA           ,times=ncompvar)
                                 , worst.rmse       = rep(-Inf         ,times=ncompvar)
                                 , worst.iata       = rep(NA_character_,times=ncompvar)
                                 , worst.n          = rep(NA           ,times=ncompvar)
                                 , stringsAsFactors = FALSE
                                 )#end data.frame
   rownames(template.table) = compvar.key
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Select all hours, and all year.                                                   #
   #---------------------------------------------------------------------------------------#
   d = which(diel.key  %in% "all.hrs")
   e = which(season.key %in% "ALL"   )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Loop over all variables.                                                         #
   #---------------------------------------------------------------------------------------#
   for (s in sequence(nsimul)){
      s.table = template.table


      for (v in sequence(ncompvar)){
         #----- Copy the variable information. --------------------------------------------#
         this.vnam     = compvar[[v]]$vnam
         this.desc     = compvar[[v]]$desc
         this.unit     = compvar[[v]]$unit
         this.sun      = compvar[[v]]$sunvar
         this.soilvar  = compvar[[v]]$soilvar
         cat("   - ",this.desc,"...","\n")

         #---------------------------------------------------------------------------------#
         #     Loop over all sites.                                                        #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata        = sites$iata[p]
            longname    = sites$desc[p]
            obser       = eft[[iata]]
            obs.now     = obser[[this.vnam]]
            nlyr        = ncol(obs.now)
            ndat        = colSums(is.finite(obs.now))
            slz.use     = ifelse(ndat==0,-obser$slz,obser$slz)
            cc          = pmin(nlyr,which.min(abs(slz.use-slz.reference)))


            #------------------------------------------------------------------------------#
            #     Grab the data for this simulation.                                       #
            #------------------------------------------------------------------------------#
            comp.now     = res[[iata]]$sim[[s]][[this.vnam]]
            n.now        = comp.now$n      [cc,d,e]
            df.now       = n.now - 1
            o.sdev.now   = sqrt(comp.now$obs.moment[cc,d,e,2])
            bias.now     = comp.now$bias   [cc,d,e] / o.sdev.now
            sigma.now    = comp.now$sigma  [cc,d,e] / o.sdev.now
            ks.stat.now  = comp.now$ks.stat[cc,d,e]
            rmse.now     = comp.now$rmse   [cc,d,e] / o.sdev.now
            #------------------------------------------------------------------------------#



            #------ Aggregate data to the total. ------------------------------------------#
            if (df.now %>% 0 & o.sdev.now %>% 0){
               s.table$n      [v] = s.table$n      [v] + n.now
               s.table$df     [v] = s.table$df     [v] + df.now
               s.table$bias   [v] = s.table$bias   [v] + n.now  * bias.now
               s.table$sigma  [v] = s.table$sigma  [v] + df.now * sigma.now^2
               s.table$ks.stat[v] = s.table$ks.stat[v] + n.now  * ks.stat.now
               if (is.finite(rmse.now) & rmse.now < s.table$best.rmse[v]){
                  s.table$best.rmse[v] = rmse.now
                  s.table$best.iata[v] = iata
                  s.table$best.n   [v] = n.now
               }#end if
               if (is.finite(rmse.now) & rmse.now > s.table$worst.rmse[v]){
                  s.table$worst.rmse[v] = rmse.now
                  s.table$worst.iata[v] = iata
                  s.table$worst.n   [v] = n.now
               }#end if
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (p in sequence(nsites))
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Normalise table.                                                               #
      #------------------------------------------------------------------------------------#
      s.table$bias    = ifelse(s.table$n  > 0,        s.table$bias/s.table$n,NA)
      s.table$sigma   = ifelse(s.table$df > 0,sqrt(s.table$sigma/s.table$df),NA)
      s.table$ks.stat = ifelse(s.table$n  > 0,     s.table$ks.stat/s.table$n,NA)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Round most of the numeric columns.                                             #
      #------------------------------------------------------------------------------------#
      s.table$bias       = round(s.table$bias      ,2)
      s.table$sigma      = round(s.table$sigma     ,2)
      s.table$ks.stat    = round(s.table$ks.stat   ,2)
      s.table$best.rmse  = round(s.table$best.rmse ,2)
      s.table$worst.rmse = round(s.table$worst.rmse,2)
      #------------------------------------------------------------------------------------#


      #----- Copy data to general table. --------------------------------------------------#
      summ.table[[simul.key[s]]] = s.table
      #------------------------------------------------------------------------------------#

   }#end for (s in sequence(nsimul))
   #---------------------------------------------------------------------------------------#
}#end if (make.summ.table)
#------------------------------------------------------------------------------------------#
