#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath" # Current directory.
srcdir         = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc" # Source  directory.
outroot        = "thisoutroot"
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
myplaces       = c("thispoly")
sasmonth       = c(2,5,8,11)
outform        = "png"          # Formats for output file.  Supported formats are:
                                 #   - "X11" - for printing on screen
                                 #   - "eps" - for postscript printing
                                 #   - "png" - for PNG printing

byeold         = TRUE           # Remove old files of the given format?

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

sasfixlimits   = FALSE          # Should I use a fixed scale for size and age-structure
                                # plots? (FALSE will set a suitable scale for each year)

ncolsfc        = 200            # Target number of colours for filled contour plots.
fcgrid         = TRUE           # Should I include a grid on the filled contour plots?

ncolshov       = 200            # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE           # Should I include a grid on the Hovmoller plots?

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.05           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.32
cex.main       = 0.8             # Scale coefficient for the title

theta           = 315.                    # Azimuth for perspective projection
phi             = 30.                     # Vertical angle for perspective projection
ltheta          = -210.                   # Azimuth angle for light
shade           = 0.125                   # Shade intensity
expz            = 0.5                     # Expansion factor for Z axis
gcol            = c("lightblue","white")  # Colours for the 50's style floor
cexmin          = 0.5                     # Minimum "head" size of the lollipop
cexmax          = 3.0                     # Maximum "head" size of the lollipop

ylnudge         = 0.05                    # Nudging factor for ylimit
ptype          = "l"                  # Type of plot
ptyped         = "p"                  # Type of plot
ptypeb         = "o"                  # Type of plot

tserdist        = TRUE          # Time series of disturbance rates

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#
#----- Time series per PFT. ---------------------------------------------------------------#
ntspft   = 20
tsplot01 = list(vnam="agbpft"       ,desc="Above ground biomass"      ,unit="kgC/m2"
               ,plt=T)
tsplot02 = list(vnam="laipft"       ,desc="Leaf area index"           ,unit="m2/m2"
               ,plt=T)
tsplot03 = list(vnam="waipft"       ,desc="Wood area index"           ,unit="m2/m2"
               ,plt=T)
tsplot04 = list(vnam="taipft"       ,desc="Total area index"          ,unit="m2/m2"
               ,plt=T)
tsplot05 = list(vnam="bseedspft"    ,desc="Seed biomass"              ,unit="kgC/m2"
               ,plt=T)
tsplot06 = list(vnam="gpppft"       ,desc="Gross primary productivity",unit="kgC/m2/yr"
               ,plt=T)
tsplot07 = list(vnam="npppft"       ,desc="Net primary productivity"  ,unit="kgC/m2/yr"
               ,plt=T)
tsplot08 = list(vnam="mcopft"       ,desc="Maintenance costs"         ,unit="kgC/m2/yr"
               ,plt=T)
tsplot09 = list(vnam="cbapft"       ,desc="Carbon balance"            ,unit="kgC/m2/yr"
               ,plt=T)
tsplot10 = list(vnam="ldroppft"     ,desc="Leaf drop"                 ,unit="kgC/m2/yr"
               ,plt=F)
tsplot11 = list(vnam="balivepft"    ,desc="Biomass of active tissues" ,unit="kgC/m2"
               ,plt=F)
tsplot12 = list(vnam="bdeadpft"     ,desc="Structural biomass"        ,unit="kgC/m2"
               ,plt=F)
tsplot13 = list(vnam="bleafpft"     ,desc="Leaf biomass"              ,unit="kgC/m2"
               ,plt=F)
tsplot14 = list(vnam="brootpft"     ,desc="Root biomass"              ,unit="kgC/m2"
               ,plt=F)
tsplot15 = list(vnam="bswoodpft"    ,desc="Sapwood biomass"           ,unit="kgC/m2"
               ,plt=F)
tsplot16 = list(vnam="bstorepft"    ,desc="Storage biomass"           ,unit="kgC/m2"
               ,plt=F)
tsplot17 = list(vnam="basareapft"   ,desc="Basal area"                ,unit="m2/m2"
               ,plt=T)
tsplot18 = list(vnam="leafresppft"  ,desc="Leaf respiration"          ,unit="m2/m2"
               ,plt=T)
tsplot19 = list(vnam="rootresppft"  ,desc="Root respiration"          ,unit="m2/m2"
               ,plt=T)
tsplot20 = list(vnam="growthresppft",desc="Growth respiration"        ,unit="m2/m2"
               ,plt=T)
#----- Time series per PFT. ---------------------------------------------------------------#
ntslu    = 6
lsplot01 = list(vnam="agblu"     ,desc="Above ground biomass"      ,unit="kgC/m2"    ,plt=F)
lsplot02 = list(vnam="lailu"     ,desc="Leaf area index"           ,unit="m2/m2"     ,plt=F)
lsplot03 = list(vnam="gpplu"     ,desc="Gross primary productivity",unit="kgC/m2/yr" ,plt=F)
lsplot04 = list(vnam="npplu"     ,desc="Net primary productivity"  ,unit="kgC/m2/yr" ,plt=F)
lsplot05 = list(vnam="arealu"    ,desc="Fraction of area"          ,unit=""          ,plt=F)
lsplot06 = list(vnam="basarealu" ,desc="Basal area"                ,unit="m2/m2"     ,plt=F)
#----- Size (DBH) and age structure of cohort level variables. ----------------------------#
npsas  = 30
psas01 = list(vnam="lightco"  ,desc="Light level"               ,unit="--"          , plt=T)
psas02 = list(vnam="beamextco",desc="Downward direct light"     ,unit="--"          , plt=T)
psas03 = list(vnam="diffextco",desc="Downward diffuse light"    ,unit="--"          , plt=T)
psas04 = list(vnam="parlco"   ,desc="Absorbed PAR "             ,unit="W/m2"        , plt=T)
psas05 = list(vnam="lambdaco" ,desc="Light extinction"          ,unit="m2/m2"       , plt=T)
psas06 = list(vnam="gppco"    ,desc="Gross primary productivity",unit="kgC/plant/yr", plt=T)
psas07 = list(vnam="respco"   ,desc="Total plant respiration"   ,unit="kgC/plant/yr", plt=T)
psas08 = list(vnam="nppco"    ,desc="Net primary productivity"  ,unit="kgC/plant/yr", plt=T)
psas09 = list(vnam="cbrbarco" ,desc="Relative carbon balance"   ,unit="--"          , plt=T)
psas10 = list(vnam="cbalco"   ,desc="Carbon balance"            ,unit="kgC/plant/yr", plt=T)
psas11 = list(vnam="mcostco"  ,desc="Maintenance costs"         ,unit="kgC/plant/yr", plt=T)
psas12 = list(vnam="ncbmortco",desc="Mortality due to Neg. CB"  ,unit="1/yr"        , plt=T)
psas13 = list(vnam="agbco"    ,desc="Above-ground biomass"      ,unit="kgC/plant"   , plt=T)
psas14 = list(vnam="fsoco"    ,desc="Fraction of open stomata"  ,unit="--"          , plt=T)
psas15 = list(vnam="nplantco" ,desc="Plant density"             ,unit="plant/m2"    , plt=T)
psas16 = list(vnam="laico"    ,desc="Leaf area index"           ,unit="m2/m2"       , plt=T)
psas17 = list(vnam="waico"    ,desc="Wood area index"           ,unit="m2/m2"       , plt=T)
psas18 = list(vnam="taico"    ,desc="Total area index"          ,unit="m2/m2"       , plt=T)
psas19 = list(vnam="demandco" ,desc="Water demand"              ,unit="kg/m2_l/day" , plt=T)
psas20 = list(vnam="supplyco" ,desc="Water supply"              ,unit="kg/m2_l/day" , plt=T)
psas21 = list(vnam="heightco" ,desc="Cohort height"             ,unit="m"           , plt=T)
psas22 = list(vnam="gpplco"   ,desc="Gross primary productivity",unit="kgC/m2lf/yr" , plt=T)
psas23 = list(vnam="baco"     ,desc="Basal area"                ,unit="m2"          , plt=T)
psas24 = list(vnam="ldropco"  ,desc="Leaf drop"                 ,unit="kgC/plant/yr", plt=T)
psas25 = list(vnam="baliveco" ,desc="Biomass of active tissues" ,unit="kgC/plant"   , plt=T)
psas26 = list(vnam="bdeadco"  ,desc="Structural biomass"        ,unit="kgC/plant"   , plt=T)
psas27 = list(vnam="bleafco"  ,desc="Leaf biomass"              ,unit="kgC/plant"   , plt=T)
psas28 = list(vnam="brootco"  ,desc="Root biomass"              ,unit="kgC/plant"   , plt=T)
psas29 = list(vnam="bswoodco" ,desc="Sapwood biomass"           ,unit="kgC/plant"   , plt=T)
psas30 = list(vnam="bstoreco" ,desc="Storage biomass"           ,unit="kgC/plant"   , plt=T)
#----- Time series of some variables per PFT, with size or age distribution. --------------#
nfcpft  = 16
fcpft01 = list(vnam="agbpftdbh",desc="Above-ground biomass"      ,unit="kgC/m2"
                               ,cls="dbh",csch="iatlas",plt=F)
fcpft02 = list(vnam="laipftdbh",desc="Leaf area index"           ,unit="m2/m2"
                               ,cls="dbh",csch="iatlas",plt=F)
fcpft03 = list(vnam="waipftdbh",desc="Wood area index"           ,unit="m2/m2"
                               ,cls="dbh",csch="iatlas",plt=F)
fcpft04 = list(vnam="taipftdbh",desc="Total area index"          ,unit="m2/m2"
                               ,cls="dbh",csch="iatlas",plt=F)
fcpft05 = list(vnam="gpppftdbh",desc="Gross primary productivity",unit="kgC/m2/yr"
                               ,cls="dbh",csch="muitas",plt=F)
fcpft06 = list(vnam="npppftdbh",desc="Net primary productivity"  ,unit="kgC/m2/yr"
                               ,cls="dbh",csch="muitas",plt=F)
fcpft07 = list(vnam="mcopftdbh",desc="Maintenance costs"         ,unit="kgC/m2/yr"
                               ,cls="dbh",csch="muitas",plt=F)
fcpft08 = list(vnam="cbapftdbh",desc="Carbon balance"            ,unit="kgC/m2/yr"
                               ,cls="dbh",csch="muitas",plt=F)
fcpft09 = list(vnam="agbpftage",desc="Above-ground biomass"      ,unit="kgC/m2"
                               ,cls="age",csch="iatlas",plt=F)
fcpft10 = list(vnam="laipftage",desc="Leaf area index"           ,unit="m2/m2"
                               ,cls="age",csch="iatlas",plt=F)
fcpft11 = list(vnam="waipftage",desc="Wood area index"           ,unit="m2/m2"
                               ,cls="age",csch="iatlas",plt=F)
fcpft12 = list(vnam="taipftage",desc="Total area index"          ,unit="m2/m2"
                               ,cls="age",csch="iatlas",plt=F)
fcpft13 = list(vnam="gpppftage",desc="Gross primary productivity",unit="kgC/m2/yr"
                               ,cls="age",csch="muitas",plt=F)
fcpft14 = list(vnam="npppftage",desc="Net primary productivity"  ,unit="kgC/m2/yr"
                               ,cls="age",csch="muitas",plt=F)
fcpft15 = list(vnam="mcopftage",desc="Maintenance costs"         ,unit="kgC/m2/yr"
                               ,cls="age",csch="muitas",plt=F)
fcpft16 = list(vnam="cbapftage",desc="Carbon balance"            ,unit="kgC/m2/yr"
                               ,cls="age",csch="muitas",plt=F)
#----- Box plots --------------------------------------------------------------------------#
nbox = 24
bplot01 = list(vnam="gpp"        ,desc="Gross Primary productivity"      ,unit="kgC/m2/yr"
                                 ,plt=T)
bplot02 = list(vnam="plresp"     ,desc="Plant respiration"               ,unit="kgC/m2/yr"
                                 ,plt=T)
bplot03 = list(vnam="hetresp"    ,desc="Heterotrophic respiration"       ,unit="kgC/m2/yr"
                                 ,plt=T)
bplot04 = list(vnam="nep"        ,desc="Net ecosystem production"        ,unit="kgC/m2/yr"
                                 ,plt=T)
bplot05 = list(vnam="sens"       ,desc="Sensible heat flux"              ,unit="W/m2"
                                 ,plt=T)
bplot06 = list(vnam="evap"       ,desc="Evaporation  "                   ,unit="kg/m2/day"
                                 ,plt=T)
bplot07 = list(vnam="transp"     ,desc="Transpiration"                   ,unit="kg/m2/day"
                                 ,plt=T)
bplot08 = list(vnam="atm.temp"   ,desc="Atmospheric temperature"         ,unit="degC"
                                 ,plt=T)
bplot09 = list(vnam="atm.shv"    ,desc="Atmospheric specific humidity"   ,unit="g/kg"
                                 ,plt=T)
bplot10 = list(vnam="atm.co2"    ,desc="Atmospheric CO2 mixing ratio"    ,unit="umol/mol"
                                 ,plt=T)
bplot11 = list(vnam="can.temp"   ,desc="Canopy air temperature"          ,unit="degC"
                                 ,plt=T)
bplot12 = list(vnam="can.shv"    ,desc="Canopy air specific humidity"    ,unit="g/kg"
                                 ,plt=T)
bplot13 = list(vnam="can.co2"    ,desc="Canopy air CO2 mixing ratio"     ,unit="umol/mol"
                                 ,plt=T)
bplot14 = list(vnam="rain"       ,desc="Total monthly precipitation"     ,unit="mm"
                                 ,plt=T)
bplot15 = list(vnam="leaf.temp"  ,desc="Leaf temperature"                ,unit="degC"
                                 ,plt=T)
bplot16 = list(vnam="wood.temp"  ,desc="Wood temperature"                ,unit="degC"
                                 ,plt=T)
bplot17 = list(vnam="gnd.temp"   ,desc="Ground temperature"              ,unit="degC"
                                 ,plt=T)
bplot18 = list(vnam="et"         ,desc="Evapotranspiration"              ,unit="kg/m2/day"
                                 ,plt=T)
bplot19 = list(vnam="fs.open"    ,desc="Fraction of open stomata"        ,unit="---"
                                 ,plt=T)
bplot20 = list(vnam="rshort"     ,desc="Downward Shortwave radiation"    ,unit="W/m2"
                                 ,plt=T)
bplot21 = list(vnam="rlong"      ,desc="Downward Longwave radiation"     ,unit="W/m2"
                                 ,plt=T)
bplot22 = list(vnam="gnd.shv"    ,desc="Ground specific humidity"        ,unit="g/kg"
                                 ,plt=T)
bplot23 = list(vnam="rshort.beam",desc="Direct incident SW radiation"    ,unit="W/m2"
                                 ,plt=T)
bplot24 = list(vnam="rshort.diff",desc="Diffuse incident SW radiation"   ,unit="W/m2"
                                 ,plt=T)
#----- Similar to Hovmoller diagrams. -----------------------------------------------------#
nhov = 37
hovdi01 = list(vnam="gpp"        ,desc="Gross Primary productivity"      ,unit="kgC/m2/yr"
                                 ,csch="atlas"                           ,plt=T)
hovdi02 = list(vnam="plresp"     ,desc="Plant respiration"               ,unit="kgC/m2/yr"
                                 ,csch="muitas"                          ,plt=T)
hovdi03 = list(vnam="hetresp"    ,desc="Heterotrophic respiration"       ,unit="kgC/m2/yr"
                                 ,csch="muitas"                          ,plt=T)
hovdi04 = list(vnam="npp"        ,desc="Net primary production"          ,unit="kgC/m2/yr"
                                 ,csch="muitas"                          ,plt=T)
hovdi05 = list(vnam="sens"       ,desc="Sensible heat flux"              ,unit="W/m2"
                                 ,csch="muitas"                          ,plt=T)
hovdi06 = list(vnam="evap"       ,desc="Evaporation"                     ,unit="kg/m2/day"
                                 ,csch="imuitas"                         ,plt=T)
hovdi07 = list(vnam="transp"     ,desc="Transpiration"                   ,unit="kg/m2/day"
                                 ,csch="imuitas"                         ,plt=T)
hovdi08 = list(vnam="atm.temp"   ,desc="Atmospheric temperature"         ,unit="degC"
                                 ,csch="muitas"                          ,plt=T)
hovdi09 = list(vnam="atm.shv"    ,desc="Atmospheric specific humidity"   ,unit="g/kg"
                                 ,csch="imuitas"                         ,plt=T)
hovdi10 = list(vnam="atm.co2"    ,desc="Atmospheric CO2 mixing ratio"    ,unit="umol/mol"
                                 ,csch="muitas"                          ,plt=F)
hovdi11 = list(vnam="can.temp"   ,desc="Canopy air temperature"          ,unit="degC"
                                 ,csch="muitas"                          ,plt=T)
hovdi12 = list(vnam="can.shv"    ,desc="Canopy air specific humidity"    ,unit="g/kg"
                                 ,csch="imuitas"                         ,plt=T)
hovdi13 = list(vnam="can.co2"    ,desc="Canopy air CO2 mixing ratio"     ,unit="umol/mol"
                                 ,csch="muitas"                          ,plt=T)
hovdi14 = list(vnam="rain"       ,desc="Total monthly precipitation"     ,unit="mm"
                                 ,csch="imuitas"                         ,plt=T)
hovdi15 = list(vnam="leaf.temp"  ,desc="Leaf temperature"                ,unit="degC"
                                 ,csch="muitas"                          ,plt=T)
hovdi16 = list(vnam="wood.temp"  ,desc="Wood temperature"                ,unit="degC"
                                 ,csch="muitas"                          ,plt=T)
hovdi17 = list(vnam="gnd.temp"   ,desc="Ground temperature"              ,unit="degC"
                                 ,csch="muitas"                          ,plt=T)
hovdi18 = list(vnam="gnd.shv"    ,desc="Ground specific humidity"        ,unit="g/kg"
                                 ,csch="imuitas"                         ,plt=T)
hovdi19 = list(vnam="workload"   ,desc="Workload"                        ,unit="steps/day"
                                 ,csch="muitas"                          ,plt=T)
hovdi20 = list(vnam="et"         ,desc="Evapotranspiration"              ,unit="kg/m2/day"
                                 ,csch="imuitas"                         ,plt=T)
hovdi21 = list(vnam="fs.open"    ,desc="Fraction of open stomata"        ,unit="---"
                                 ,csch="imuitas"                         ,plt=T)
hovdi22 = list(vnam="specwork"   ,desc="Specific workload"           ,unit="steps/patch/day"
                                 ,csch="muitas"                          ,plt=T)
hovdi23 = list(vnam="wflxgc"     ,desc="Ground evaporation"              ,unit="kg/m2/day"
                                 ,csch="imuitas"                         ,plt=T)
hovdi24 = list(vnam="wflxlc"     ,desc="Leaf evaporation"                ,unit="kg/m2/day"
                                 ,csch="imuitas"                         ,plt=T)
hovdi25 = list(vnam="wflxwc"     ,desc="Wood evaporation"                ,unit="kg/m2/day"
                                 ,csch="imuitas"                         ,plt=T)
hovdi26 = list(vnam="nep"        ,desc="Net ecosystem production"        ,unit="kgC/m2/yr"
                                 ,csch="muitas"                          ,plt=T)
hovdi27 = list(vnam="nee"        ,desc="Net ecosystem exchange"          ,unit="kgC/m2/yr"
                                 ,csch="imuitas"                         ,plt=T)
hovdi28 = list(vnam="cba"         ,desc="Carbon balance"                 ,unit="kgC/m2/day"
                                 ,csch="atlas"                           ,plt=T)
hovdi29 = list(vnam="mco"        ,desc="Maintenance costs"               ,unit="kgC/m2/yr"
                                 ,csch="iatlas"                          ,plt=T)
hovdi30 = list(vnam="ldrop"      ,desc="Leaf drop"                       ,unit="kgC/m2/yr"
                                 ,csch="iatlas"                          ,plt=T)
hovdi31 = list(vnam="rshort"     ,desc="Downward shortwave radiation"    ,unit="W/m2"
                                 ,csch="icloudy"                         ,plt=T)
hovdi32 = list(vnam="rlong"      ,desc="Downward longwave radiation"     ,unit="W/m2"
                                 ,csch="cloudy"                          ,plt=T)
hovdi33 = list(vnam="rshort.gnd" ,desc="Abs. gnd. shortwave radiation"   ,unit="W/m2"
                                 ,csch="icloudy"                         ,plt=T)
hovdi34 = list(vnam="rlong.gnd"  ,desc="Abs. gnd. longwave radiation"    ,unit="W/m2"
                                 ,csch="cloudy"                          ,plt=T)
hovdi35 = list(vnam="rshort.beam",desc="Direct incident SW radiation"    ,unit="W/m2"
                                 ,csch="icloudy"                         ,plt=T)
hovdi36 = list(vnam="rshort.diff",desc="Diffuse incident SW radiation"   ,unit="W/m2"
                                 ,csch="icloudy"                         ,plt=T)
hovdi37 = list(vnam="albedo"     ,desc="SW albedo"                       ,unit="---"
                                 ,csch="muitas"                          ,plt=T)
#----- Time series with several variables in it. ------------------------------------------#
ntser=10
tser01 = list(vnam   = c("gpp","plresp","hetresp","nee")
             ,desc   = c("GPP","Plant resp.","Het. resp.","NEE")
             ,colour = c("forestgreen","chartreuse","sienna","deepskyblue")
             ,lwd    = c(2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "carbflux"
             ,theme  = "Carbon fluxes"
             ,unit   = "kgC/m2/yr"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser02 = list(vnam   = c("rshort","rlong","rshort.gnd","latent","sens")
             ,desc   = c("Down SW","Down LW","Abs. Grnd","Latent","Sensible")
             ,colour = c("goldenrod","lawngreen","purple4","steelblue","firebrick")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "eneflux"
             ,theme  = "Energy fluxes"
             ,unit   = "W/m2"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser03 = list(vnam   = c("wflxgc","et","wflxlc","wflxwc","transp")
             ,desc   = c("Ground->Canopy","Canopy->Atm","Leaf->Canopy"
                        ,"Wood->Canopy","Transpiration")
             ,colour = c("firebrick","midnightblue","chartreuse"
                        ,"goldenrod","darkolivegreen")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "h2oflux"
             ,theme  = "Water fluxes"
             ,unit   = "kg/m2/day"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser04 = list(vnam   = c("atm.temp","can.temp","leaf.temp","wood.temp","gnd.temp")
             ,desc   = c("Atmosphere","Canopy air","Leaf","Wood","Ground")
             ,colour = c("deepskyblue","gray21","chartreuse","goldenrod","sienna")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "temperature"
             ,theme  = "Temperature"
             ,unit   = "degC"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser05 = list(vnam   = c("atm.shv","can.shv","gnd.shv")
             ,desc   = c("Atmosphere","Canopy air","Ground")
             ,colour = c("deepskyblue","gray21","sienna")
             ,lwd    = c(2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "h2ovapour"
             ,theme  = "Water vapour mixing ratio"
             ,unit   = "g/kg"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser06 = list(vnam   = c("atm.co2","can.co2")
             ,desc   = c("Atmosphere","Canopy air")
             ,colour = c("deepskyblue","lawngreen")
             ,lwd    = c(2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "co2"
             ,theme  = "CO2 mixing ratio"
             ,unit   = "umol/mol"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser07 = list(vnam   = c("rain")
             ,desc   = c("Precipitation")
             ,colour = c("midnightblue")
             ,lwd    = c(2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "prec"
             ,theme  = "Precipitation rate"
             ,unit   = "mm/hr"
             ,legpos = "bottomleft"
             ,plt    = TRUE)
tser08 = list(vnam   = c("npat.global")
             ,desc   = c("Patch count")
             ,colour = c("darkorange3")
             ,lwd    = c(2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "npatch"
             ,theme  = "Total number of patches"
             ,unit   = "---"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser09 = list(vnam   = c("ncoh.global")
             ,desc   = c("Cohort count")
             ,colour = c("limegreen")
             ,lwd    = c(2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "ncohort"
             ,theme  = "Total number of cohorts"
             ,unit   = "---"
             ,legpos = "topleft"
             ,plt    = TRUE)
tser10 = list(vnam   = c("plresp","hetresp","leaf.resp","root.resp","growth.resp")
             ,desc   = c("Plant resp.","Het. resp.","Leaf resp.","Root resp."
                        ,"Growth resp.")
             ,colour = c("midnightblue","steelblue","forestgreen","darkorange3"
                        ,"goldenrod")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "respiration"
             ,theme  = "Respiration fluxes"
             ,unit   = "kgC/m2/yr"
             ,legpos = "topleft"
             ,plt    = TRUE)
#------------------------------------------------------------------------------------------#




#----- "Climatology of the mean diurnal cycle with several variables in it. ---------------#
nclim=12
clim01 = list(vnam   = c("gpp","plresp","hetresp","nep","nee")
             ,desc   = c("GPP","Plant resp.","Het. resp.","NEP","NEE")
             ,colour = c("chartreuse","goldenrod","sienna","forestgreen","deepskyblue")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "carbflux"
             ,theme  = "Carbon fluxes"
             ,unit   = "kgC/m2/yr"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim02 = list(vnam   = c("rshort","rlong","rlongup","latent","sens")
             ,desc   = c("Down SW","Down LW","Up LW","Latent","Sensible")
             ,colour = c("goldenrod","lawngreen","purple4","midnightblue","firebrick")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "eneflux"
             ,theme  = "Energy fluxes"
             ,unit   = "W/m2"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim03 = list(vnam   = c("wflxgc","et","wflxlc","wflxwc","transp")
             ,desc   = c("Ground->Canopy","Canopy->Air","Leaf->Canopy"
                        ,"Wood->Canopy","Transpiration")
             ,colour = c("firebrick","midnightblue","chartreuse"
                        ,"goldenrod","darkolivegreen")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "h2oflux"
             ,theme  = "Water fluxes"
             ,unit   = "kg/m2/day"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim04 = list(vnam   = c("atm.temp","can.temp","leaf.temp","wood.temp","gnd.temp")
             ,desc   = c("Atmosphere","Canopy air","Leaf","Wood","Ground")
             ,colour = c("deepskyblue","gray21","chartreuse","goldenrod","sienna")
             ,lwd    = c(2.5,2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "temperature"
             ,theme  = "Temperature"
             ,unit   = "degC"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim05 = list(vnam   = c("atm.shv","can.shv","gnd.shv")
             ,desc   = c("Atmosphere","Canopy air","Ground")
             ,colour = c("deepskyblue","gray21","sienna")
             ,lwd    = c(2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "h2ovapour"
             ,theme  = "Water vapour mixing ratio"
             ,unit   = "g/kg"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim06 = list(vnam   = c("atm.co2","can.co2")
             ,desc   = c("Atmosphere","Canopy air")
             ,colour = c("deepskyblue","lawngreen")
             ,lwd    = c(2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "co2"
             ,theme  = "CO2 mixing ratio"
             ,unit   = "umol/mol"
             ,legpos = "bottomleft"
             ,plt    = TRUE)
clim07 = list(vnam   = c("rain")
             ,desc   = c("Precipitation")
             ,colour = c("midnightblue")
             ,lwd    = c(2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "prec"
             ,theme  = "Precipitation rate"
             ,unit   = "mm/hr"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim08 = list(vnam   = c("fs.open")
             ,desc   = c("Fraction of open stomata")
             ,colour = c("green3")
             ,lwd    = c(2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "fsopen"
             ,theme  = "Fraction of open stomata"
             ,unit   = "---"
             ,legpos = "bottomleft"
             ,plt    = TRUE)
clim09 = list(vnam   = c("rshort","rshort.beam","rshort.diff","rshort.gnd")
             ,desc   = c("Down Top canopy","Beam","Diffuse","Abs. Ground")
             ,colour = c("deepskyblue","goldenrod","gray33","firebrick")
             ,lwd    = c(2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "rshort"
             ,theme  = "Short wave radiation"
             ,unit   = "W/m2"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim10 = list(vnam   = c("rlong","rlongup","rlong.gnd")
             ,desc   = c("Down Top canopy","Upward LW","Abs. Ground")
             ,colour = c("deepskyblue","firebrick","goldenrod")
             ,lwd    = c(2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "rlong"
             ,theme  = "Long wave radiation"
             ,unit   = "W/m2"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim11 = list(vnam   = c("albedo","albedo.beam","albedo.diff","rlong.albedo")
             ,desc   = c("SW Albedo (Net)","SW Albedo (Beam)","SW Albedo (Diff)"
                        ,"LW Albedo")
             ,colour = c("deepskyblue","goldenrod","gray33","firebrick")
             ,lwd    = c(2.5,2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "albedo"
             ,theme  = "Albedo"
             ,unit   = "---"
             ,legpos = "topleft"
             ,plt    = TRUE)
clim12 = list(vnam   = c("hflxgc","sens","hflxlc","hflxwc")
             ,desc   = c("Ground->Canopy","Canopy->Air","Leaf->Canopy","Wood->Canopy")
             ,colour = c("firebrick","midnightblue","chartreuse","goldenrod")
             ,lwd    = c(2.5,2.5,2.5)
             ,type   = ptype
             ,plog   = ""
             ,prefix = "heatflux"
             ,theme  = "Sensible heat fluxes"
             ,unit   = "kg/m2/day"
             ,legpos = "topleft"
             ,plt    = TRUE)
#------------------------------------------------------------------------------------------#


#----- Similar to Hovmoller diagrams, but with month/year in the y axis and hour in x. ----#
nhdcyc  = 29
hdcyc01 = list(vnam="gpp"       ,desc="Gross Primary productivity"      ,unit="kgC/m2/yr"
                                ,csch="atlas"                           ,plt=T)
hdcyc02 = list(vnam="plresp"    ,desc="Plant respiration"               ,unit="kgC/m2/yr"
                                ,csch="muitas"                          ,plt=T)
hdcyc03 = list(vnam="hetresp"   ,desc="Heterotrophic respiration"       ,unit="kgC/m2/yr"
                                ,csch="muitas"                          ,plt=T)
hdcyc04 = list(vnam="sens"      ,desc="Sensible heat flux"              ,unit="W/m2"     
                                ,csch="muitas"                          ,plt=T)
hdcyc05 = list(vnam="evap"      ,desc="Evaporation"                     ,unit="kg/m2/day"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc06 = list(vnam="transp"    ,desc="Transpiration"                   ,unit="kg/m2/day"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc07 = list(vnam="atm.temp"  ,desc="Atmospheric temperature"         ,unit="degC"     
                                ,csch="muitas"                          ,plt=T)
hdcyc08 = list(vnam="atm.shv"   ,desc="Atmospheric specific humidity"   ,unit="g/kg"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc09 = list(vnam="atm.co2"   ,desc="Atmospheric CO2 mixing ratio"    ,unit="umol/mol" 
                                ,csch="muitas"                          ,plt=F)
hdcyc10 = list(vnam="can.temp"  ,desc="Canopy air temperature"          ,unit="degC"     
                                ,csch="muitas"                          ,plt=T)
hdcyc11 = list(vnam="can.shv"   ,desc="Canopy air specific humidity"    ,unit="g/kg"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc12 = list(vnam="can.co2"   ,desc="Canopy air CO2 mixing ratio"     ,unit="umol/mol" 
                                ,csch="muitas"                          ,plt=T)
hdcyc13 = list(vnam="rain"      ,desc="Total monthly precipitation"     ,unit="mm"       
                                ,csch="imuitas"                         ,plt=T)
hdcyc14 = list(vnam="leaf.temp" ,desc="Leaf temperature"                ,unit="degC"     
                                ,csch="muitas"                          ,plt=T)
hdcyc15 = list(vnam="wood.temp" ,desc="Wood temperature"                ,unit="degC"     
                                ,csch="muitas"                          ,plt=T)
hdcyc16 = list(vnam="gnd.temp"  ,desc="Ground temperature"              ,unit="degC"     
                                ,csch="muitas"                          ,plt=T)
hdcyc17 = list(vnam="gnd.shv"   ,desc="Ground specific humidity"        ,unit="g/kg"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc18 = list(vnam="et"        ,desc="'Evapotranspiration"             ,unit="kg/m2/day"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc19 = list(vnam="fs.open"   ,desc="Fraction of open stomata"        ,unit="---"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc20 = list(vnam="wflxgc"    ,desc="Ground evaporation"              ,unit="kg/m2/day"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc21 = list(vnam="wflxlc"    ,desc="Leaf evaporation"                ,unit="kg/m2/day"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc22 = list(vnam="wflxwc"    ,desc="Wood evaporation"                ,unit="kg/m2/day"     
                                ,csch="imuitas"                         ,plt=T)
hdcyc23 = list(vnam="nep"       ,desc="Net ecosystem production"        ,unit="kgC/m2/yr"
                                ,csch="muitas"                          ,plt=T)
hdcyc24 = list(vnam="nee"       ,desc="Net ecosystem exchange"          ,unit="kgC/m2/yr"
                                ,csch="imuitas"                         ,plt=T)
hdcyc25 = list(vnam="rshort"    ,desc="Downward shortwave radiation"    ,unit="W/m2"
                                ,csch="icloudy"                         ,plt=T)
hdcyc26 = list(vnam="rlong"     ,desc="Downward longwave radiation"     ,unit="W/m2"
                                ,csch="cloudy"                          ,plt=T)
hdcyc27 = list(vnam="rshort.gnd",desc="Abs. gnd. shortwave radiation"   ,unit="W/m2"
                                ,csch="icloudy"                         ,plt=T)
hdcyc28 = list(vnam="rlong.gnd" ,desc="Abs. gnd. longwave radiation"    ,unit="W/m2"
                                ,csch="cloudy"                          ,plt=T)
hdcyc29 = list(vnam="albedo"    ,desc="Shortwave albedo"                ,unit="---"
                                ,csch="muitas"                          ,plt=T)
#------------------------------------------------------------------------------------------#




#----- Comparison between observations and model averages. --------------------------------#
ncompdcyc = 12
compdcyc01 = list( vnam   = "gpp"
                 , desc   = "Gross Primary Productivity"
                 , unit   = "kgC/m2/yr"
                 , plotsd = TRUE
                 , colour = c("darkgreen","gray14")
                 , errcol = c("chartreuse","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc02 = list( vnam   = "nee"
                 , desc   = "Net Ecosystem Exchange"
                 , unit   = "kgC/m2/yr"
                 , plotsd = TRUE
                 , colour = c("darkgreen","gray14")
                 , errcol = c("chartreuse","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "bottomleft"
                 , plt    = TRUE)
compdcyc03 = list( vnam   = "et"
                 , desc   = "Evapotranspiration"
                 , unit   = "kgH2O/m2/day"
                 , plotsd = TRUE
                 , colour = c("midnightblue","gray14")
                 , errcol = c("steelblue","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc04 = list( vnam   = "sens"
                 , desc   = "Sensible heat flux"
                 , unit   = "W/m2"
                 , plotsd = TRUE
                 , colour = c("darkorange3","gray14")
                 , errcol = c("gold","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc05 = list( vnam   = "rain"
                 , desc   = "Precipitation rate"
                 , unit   = "kgH2O/m2/day"
                 , plotsd = FALSE
                 , colour = c("midnightblue","gray14")
                 , errcol = c("steelblue","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc06 = list( vnam   = "atm.temp"
                 , desc   = "Air temperature"
                 , unit   = "C"
                 , plotsd = TRUE
                 , colour = c("darkorange3","gray14")
                 , errcol = c("gold","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc07 = list( vnam   = "rshort"
                 , desc   = "Shortwave radiation"
                 , unit   = "W/m2"
                 , plotsd = TRUE
                 , colour = c("darkorange3","gray14")
                 , errcol = c("gold","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc08 = list( vnam   = "rlong"
                 , desc   = "Longwave radiation"
                 , unit   = "W/m2"
                 , plotsd = TRUE
                 , colour = c("midnightblue","gray14")
                 , errcol = c("steelblue","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc09 = list( vnam   = "atm.shv"
                 , desc   = "Air specific humidity"
                 , unit   = "g/kg"
                 , plotsd = TRUE
                 , colour = c("midnightblue","gray14")
                 , errcol = c("steelblue","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc10 = list( vnam   = "atm.co2"
                 , desc   = "Air CO2 mixing ratio"
                 , unit   = "g/kg"
                 , plotsd = TRUE
                 , colour = c("darkgreen","gray14")
                 , errcol = c("chartreuse","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc11 = list( vnam   = "atm.prss"
                 , desc   = "Air pressure"
                 , unit   = "hPa"
                 , plotsd = TRUE
                 , colour = c("purple4","gray14")
                 , errcol = c("orchid","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
compdcyc12 = list( vnam   = "atm.vels"
                 , desc   = "Wind speed"
                 , unit   = "m/s"
                 , plotsd = TRUE
                 , colour = c("purple4","gray14")
                 , errcol = c("orchid","gray61")
                 , angle  = c(45,-45)
                 , dens   = c(40, 40)
                 , lwd    = c(2.5,2.5)
                 , shwd   = c(1.0,1.0)
                 , type   = "o"
                 , plog   = ""
                 , legpos = "topleft"
                 , plt    = TRUE)
#------------------------------------------------------------------------------------------#




#----- Comparison between observations and model averages. --------------------------------#
ncompmmean = 12
compmmean01 = list( vnam   = "gpp"
                  , desc   = "Gross Primary Productivity"
                  , unit   = "kgC/m2/yr"
                  , plotsd = TRUE
                  , colour = c("darkgreen","gray14")
                  , errcol = c("chartreuse","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean02 = list( vnam   = "nee"
                  , desc   = "Net Ecosystem Exchange"
                  , unit   = "kgC/m2/yr"
                  , plotsd = TRUE
                  , colour = c("darkgreen","gray14")
                  , errcol = c("chartreuse","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean03 = list( vnam   = "et"
                  , desc   = "Evapotranspiration"
                  , unit   = "kgH2O/m2/day"
                  , plotsd = TRUE
                  , colour = c("midnightblue","gray14")
                  , errcol = c("steelblue","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean04 = list( vnam   = "sens"
                  , desc   = "Sensible heat"
                  , unit   = "W/m2"
                  , plotsd = TRUE
                  , colour = c("darkorange3","gray14")
                  , errcol = c("gold","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean05 = list( vnam   = "rain"
                  , desc   = "Precipitation rate"
                  , unit   = "kgH2O/m2/day"
                  , plotsd = FALSE
                  , colour = c("midnightblue","gray14")
                  , errcol = c("steelblue","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean06 = list( vnam   = "atm.temp"
                  , desc   = "Air temperature"
                  , unit   = "C"
                  , plotsd = TRUE
                  , colour = c("darkorange3","gray14")
                  , errcol = c("gold","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean07 = list( vnam   = "rshort"
                  , desc   = "Shortwave radiation"
                  , unit   = "W/m2"
                  , plotsd = TRUE
                  , colour = c("darkorange3","gray14")
                  , errcol = c("gold","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean08 = list( vnam   = "rlong"
                  , desc   = "Longwave radiation"
                  , unit   = "W/m2"
                  , plotsd = TRUE
                  , colour = c("midnightblue","gray14")
                  , errcol = c("steelblue","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean09 = list( vnam   = "atm.shv"
                  , desc   = "Air specific humidity"
                  , unit   = "g/kg"
                  , plotsd = TRUE
                  , colour = c("midnightblue","gray14")
                  , errcol = c("steelblue","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean10 = list( vnam   = "atm.co2"
                  , desc   = "Air CO2 mixing ratio"
                  , unit   = "g/kg"
                  , plotsd = TRUE
                  , colour = c("darkgreen","gray14")
                  , errcol = c("chartreuse","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean11 = list( vnam   = "atm.prss"
                  , desc   = "Air pressure"
                  , unit   = "hPa"
                  , plotsd = TRUE
                  , colour = c("purple4","gray14")
                  , errcol = c("orchid","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
compmmean12 = list( vnam   = "atm.vels"
                  , desc   = "Wind speed"
                  , unit   = "m/s"
                  , plotsd = TRUE
                  , colour = c("purple4","gray14")
                  , errcol = c("orchid","gray61")
                  , angle  = c(45,-45)
                  , dens   = c(40, 40)
                  , lwd    = c(2.5,2.5)
                  , shwd   = c(1.0,1.0)
                  , type   = "o"
                  , plog   = ""
                  , legpos = "topleft"
                  , plt    = TRUE)
#------------------------------------------------------------------------------------------#



#----- Annual mean. -----------------------------------------------------------------------#
nsoilclim  = 3
soilclim01 = list(vnam="soil.water" ,desc = "Soil moisture"
                 ,unit = "m3/m3"    ,csch = "imuitas"
                 ,pnlog=FALSE       ,plt  = TRUE)
soilclim02 = list(vnam="soil.temp"  ,desc = "Soil temperature"
                 ,unit = "C"        ,csch = "muitas"
                 ,pnlog=FALSE       ,plt  = TRUE)
soilclim03 = list(vnam="soil.mstpot",desc = "(Negative) Soil moisture potential"
                 ,unit = "m"        ,csch = "muitas"
                 ,pnlog=TRUE        ,plt  = TRUE)
#------------------------------------------------------------------------------------------#




#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(Hmisc)

#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()

#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)

#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)


#----- Load some files with functions. ----------------------------------------------------#
source(paste(srcdir,"atlas.r"      ,sep="/"))
source(paste(srcdir,"charutils.r"  ,sep="/"))
source(paste(srcdir,"census.r"     ,sep="/"))
source(paste(srcdir,"cloudy.r"     ,sep="/"))
source(paste(srcdir,"error.bar.r"  ,sep="/"))
source(paste(srcdir,"globdims.r"   ,sep="/"))
source(paste(srcdir,"locations.r"  ,sep="/"))
source(paste(srcdir,"muitas.r"     ,sep="/"))
source(paste(srcdir,"plotsize.r"   ,sep="/"))
source(paste(srcdir,"pretty.log.r" ,sep="/"))
source(paste(srcdir,"pretty.time.r",sep="/"))
source(paste(srcdir,"qapply.r"     ,sep="/"))
source(paste(srcdir,"rconstants.r" ,sep="/"))
source(paste(srcdir,"soilutils.r"  ,sep="/"))
source(paste(srcdir,"sombreado.r"  ,sep="/"))
source(paste(srcdir,"southammap.r" ,sep="/"))
source(paste(srcdir,"timeutils.r"  ,sep="/"))


#----- Load observations. -----------------------------------------------------------------#
obsrfile = paste(srcdir,"LBA_MIP.v3.RData",sep="/")
load(file=obsrfile)


#----- Define some default legend colours and names. --------------------------------------#
pftnames = c("C4 Grass","Early Tropical","Mid Tropical","Late Tropical","Temp. C3 Grass"
             ,"North Pine","South Pine","Late Conifer","Early Temperate","Mid Temperate"
             ,"Late Temperate","C3 Pasture","C3 Crop","C4 Pasture","C4 Crop"
             ,"C3 Grass","Araucaria","Total")
pftcols  = c("gold","chartreuse","limegreen","darkgreen","purple3"
            ,"deepskyblue","aquamarine","midnightblue","darkorange3","sienna"
            ,"firebrick","orchid","coral","gray45","olivedrab"
            ,"goldenrod","steelblue","gray22")

lunames = c("Agricultural","Secondary","Primary","Total")
lucols  = c("goldenrod","chartreuse","darkgreen","gray22")

distnames = c("Agr->Agr" ,"2nd->Agr" ,"Prim->Agr"
              ,"Agr->2nd" ,"2nd->2nd" ,"Prim->2nd"
              ,"Agr->Prim","2nd->Prim","Prim->Prim")
distcols  = c("gold","darkorange2","firebrick"
              ,"lightskyblue","turquoise","steelblue"
              ,"palegreen","chartreuse","forestgreen")

#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)

#------------------------------------------------------------------------------------------#
#     Using brute force to concatenate the plotting lists.                                 #
#------------------------------------------------------------------------------------------#
#----- Time series by PFT plot. -----------------------------------------------------------#
tspft = list()
for (s in 1:ntspft){
  sss           = substring(100+s,2,3)
  tsts          = paste("tsplot",sss,sep="")
  tspft[[s]]    = get(tsts)
} #end for
#----- Time series by LU plot. ------------------------------------------------------------#
tslu = list()
for (s in 1:ntslu){
  sss           = substring(100+s,2,3)
  lulu          = paste("lsplot",sss,sep="")
  tslu[[s]]     = get(lulu)
} #end for
#----- Size (DBH) and age structure cohort level variables. -------------------------------#
sasplot = list()
for (s in 1:npsas){
  sss           = substring(100+s,2,3)
  psps          = paste("psas",sss,sep="")
  sasplot[[s]]  = get(psps)
} #end for
#----- Time series of some variables per PFT, with size or age distribution. ---------------#
fcpft = list()
for (s in 1:nfcpft){
  sss         = substring(100+s,2,3)
  fcfc        = paste("fcpft",sss,sep="")
  fcpft[[s]]  = get(fcfc)
} #end for
#----- Box plots --------------------------------------------------------------------------#
bplot = list()
for (s in 1:nbox){
  sss         = substring(100+s,2,3)
  bpbp        = paste("bplot",sss,sep="")
  bplot[[s]]  = get(bpbp)
} #end for
#----- Hovmoller diagram ------------------------------------------------------------------#
hovdi = list()
for (s in 1:nhov){
  sss         = substring(100+s,2,3)
  hdhd        = paste("hovdi",sss,sep="")
  hovdi[[s]]  = get(hdhd)
} #end for
#----- Time series ------------------------------------------------------------------------#
tser      = list()
for (s in 1:ntser){
  sss         = substring(100+s,2,3)
  tsts        = paste("tser",sss,sep="")
  tser[[s]]   = get(tsts)
} #end for
#----- Time series of "mean diurnal cycle climatology" ------------------------------------#
clim      = list()
for (s in 1:nclim){
  sss        = substring(100+s,2,3)
  clcl       = paste("clim",sss,sep="")
  clim[[s]]  = get(clcl)
} #end for
#----- Hovmoller diagram ------------------------------------------------------------------#
hdcyc      = list()
for (s in 1:nhdcyc){
  sss        = substring(100+s,2,3)
  hdhd       = paste("hdcyc",sss,sep="")
  hdcyc[[s]] = get(hdhd)
} #end for
#----- Comparison of mean diurnal cycle ---------------------------------------------------#
compdcyc      = list()
for (o in 1:ncompdcyc){
  ccc          = substring(100+o,2,3)
  cpcp         = paste("compdcyc",ccc,sep="")
  compdcyc[[o]] = get(cpcp)
} #end for
#----- Comparison of monthly means --------------------------------------------------------#
compmmean      = list()
for (o in 1:ncompmmean){
  ccc            = substring(100+o,2,3)
  cmcm           = paste("compmmean",ccc,sep="")
  compmmean[[o]] = get(cmcm)
} #end for
#----- Soil profile climatology. ----------------------------------------------------------#
soilclim      = list()
for (o in 1:nsoilclim){
  ccc           = substring(100+o,2,3)
  scsc          = paste("soilclim",ccc,sep="")
  soilclim[[o]] = get(scsc)
} #end for
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi = locations(where=place,here=here,yearbeg=yearbeg,yearend=yearend
                      ,monthbeg=monthbeg)
   inpref  = thispoi$pathin
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"monthly",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz

   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outpref)) dir.create(outpref)

   #----- Print a banner to entretain the user. -------------------------------------------#
   print(paste(" + Post-processing output from ",lieu,"...",sep=""))


   #---------------------------------------------------------------------------------------#
   #     Flush all variables that will hold the data.                                      #
   #---------------------------------------------------------------------------------------#
   totmon      = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   #----- Size (DBH) and age arrays. ------------------------------------------------------#
   agbpftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   laipftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   waipftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   taipftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   gpppftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   npppftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   mcopftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   cbapftdbh      = array(data=0.,dim=c(totmon,ndbh,npft))
   agbpftage      = array(data=0.,dim=c(totmon,nage,npft))
   laipftage      = array(data=0.,dim=c(totmon,nage,npft))
   waipftage      = array(data=0.,dim=c(totmon,nage,npft))
   taipftage      = array(data=0.,dim=c(totmon,nage,npft))
   gpppftage      = array(data=0.,dim=c(totmon,nage,npft))
   npppftage      = array(data=0.,dim=c(totmon,nage,npft))
   mcopftage      = array(data=0.,dim=c(totmon,nage,npft))
   cbapftage      = array(data=0.,dim=c(totmon,nage,npft))
   #----- PFT arrays.   The "+1" column contains the total. -------------------------------#
   agbpft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   bseedspft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   laipft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   waipft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   taipft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   gpppft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   npppft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   mcopft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   cbapft         = matrix(data=0,nrow=totmon,ncol=npft+1)
   ldroppft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   balivepft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   bdeadpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   bleafpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   brootpft       = matrix(data=0,nrow=totmon,ncol=npft+1)
   bswoodpft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   bstorepft      = matrix(data=0,nrow=totmon,ncol=npft+1)
   basareapft     = matrix(data=0,nrow=totmon,ncol=npft+1)
   leafresppft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   rootresppft    = matrix(data=0,nrow=totmon,ncol=npft+1)
   growthresppft  = matrix(data=0,nrow=totmon,ncol=npft+1)

   #----- LU arrays.   The "+1" column contains the total. --------------------------------#
   agblu          = matrix(data=0,nrow=totmon,ncol=nlu+1)
   lailu          = matrix(data=0,nrow=totmon,ncol=nlu+1)
   gpplu          = matrix(data=0,nrow=totmon,ncol=nlu+1)
   npplu          = matrix(data=0,nrow=totmon,ncol=nlu+1)
   arealu         = matrix(data=0,nrow=totmon,ncol=nlu+1)
   basarealu      = matrix(data=0,nrow=totmon,ncol=nlu+1)
   #----- Miscellaneous arrays. -----------------------------------------------------------#
   dist           = array(NA,dim=c(totmon,nlu,nlu))
   #----- Polygon level vectors. ----------------------------------------------------------#
   gpp             = NULL
   plresp          = NULL
   leaf.resp       = NULL
   root.resp       = NULL
   growth.resp     = NULL
   hetresp         = NULL
   mco             = NULL
   npp             = NULL
   cba             = NULL
   ldrop           = NULL
   npp             = NULL
   nep             = NULL
   nee             = NULL
   evap            = NULL
   transp          = NULL
   atm.vels        = NULL
   atm.prss        = NULL
   atm.temp        = NULL
   can.prss        = NULL
   can.temp        = NULL
   atm.co2         = NULL
   can.co2         = NULL
   leaf.temp       = NULL
   wood.temp       = NULL
   atm.shv         = NULL
   can.shv         = NULL
   can.co2         = NULL
   sens            = NULL
   latent          = NULL
   et              = NULL
   agb             = NULL
   lai             = NULL
   wai             = NULL
   tai             = NULL
   area            = NULL
   rain            = NULL
   gnd.temp        = NULL
   gnd.shv         = NULL
   workload        = NULL
   specwork        = NULL
   fs.open         = NULL
   hflxgc          = NULL
   hflxlc          = NULL
   hflxwc          = NULL
   wflxgc          = NULL
   wflxlc          = NULL
   wflxwc          = NULL
   et              = NULL
   rshort          = NULL
   rshort.beam     = NULL
   rshort.diff     = NULL
   rlong           = NULL
   rshort.gnd      = NULL
   rlong.gnd       = NULL
   rlongup         = NULL
   albedo          = NULL
   albedo.beam     = NULL
   albedo.diff     = NULL
   rlong.albedo    = NULL
   npat.global     = NULL
   ncoh.global     = NULL
   mmsqu.gpp       = NULL
   mmsqu.plresp    = NULL
   mmsqu.leaf.resp = NULL
   mmsqu.root.resp = NULL
   mmsqu.plresp    = NULL
   mmsqu.hetresp   = NULL
   mmsqu.nee       = NULL
   mmsqu.sens      = NULL
   mmsqu.hflxlc    = NULL
   mmsqu.hflxwc    = NULL
   mmsqu.hflxgc    = NULL
   mmsqu.et        = NULL
   mmsqu.latent    = NULL
   mmsqu.wflxlc    = NULL
   mmsqu.wflxwc    = NULL
   mmsqu.wflxgc    = NULL
   mmsqu.evap      = NULL
   mmsqu.transp    = NULL

   #----- Cohort level lists. -------------------------------------------------------------#
   lightco      = list()
   beamextco    = list()
   diffextco    = list()
   parlco       = list()
   lambdaco     = list()
   gppco        = list()
   gpplco       = list()
   respco       = list()
   nppco        = list()
   cbrbarco     = list()
   cbalco       = list()
   mcostco      = list()
   ldropco      = list()
   ncbmortco    = list()
   agbco        = list()
   fsoco        = list()
   nplantco     = list()
   pftco        = list()
   dbhco        = list()
   laico        = list()
   waico        = list()
   taico        = list()
   ageco        = list()
   areaco       = list()
   demandco     = list()
   supplyco     = list()
   heightco     = list()
   baco         = list()
   baliveco     = list()
   bdeadco      = list()
   bleafco      = list()
   brootco      = list()
   bswoodco     = list()
   bstoreco     = list()

   n            = 0
   m            = 0
   thismonth    = NULL
   monnum       = NULL
   myear        = NULL

   first.time   = TRUE

   #----- Loop over years. ----------------------------------------------------------------#
   for (year in yeara:yearz){
       if (year == yeara){
          firstmonth = monthbeg
       }else{
          firstmonth = 1
       }#end if
       if (year == yearz){
          lastmonth = meszz
       }else{
          lastmonth = 12
       }#end if
       print (paste("    - Reading data from year ",year,"...",sep=""))

       #----- Loop over months. -----------------------------------------------------------#
       for (month in firstmonth:lastmonth){
          m = m + 1

          #----- Build the month and year vector. -----------------------------------------#
          monnum = c(monnum,month)
          myear  = c(myear,year)

          #----- Build the file name. -----------------------------------------------------#
          cmonth = substring(100+month,2,3)
          ddd    = daymax(month,year)
          myfile = paste(inpref,"-Q-",year,"-",cmonth,"-00-000000-g01.h5",sep="")
          #--------------------------------------------------------------------------------#



          #----- Read data and close connection immediately after. ------------------------#
          mymont = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
          #--------------------------------------------------------------------------------#


          #----- Build the time. ----------------------------------------------------------#
          thismonth = c(thismonth,chron(paste(month,1,year,sep="/"),times="0:0:0"))
          #--------------------------------------------------------------------------------#


          #----- Define the number of soil layers. ----------------------------------------#
          nzg      = mymont$NZG
          nzs      = mymont$NZS
          ndcycle  = mymont$NDCYCLE
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #      If this is the first time, allocate and initialise the mean diurnal cycle #
          # arrays.                                                                        #
          #--------------------------------------------------------------------------------#
          if (first.time){
             first.time          = FALSE

             #----- Find which soil are we solving, and save properties into soil.prop. ---#
             isoilflg   = mymont$ISOILFLG
             slz        = mymont$SLZ
             slxsand    = mymont$SLXSAND
             slxclay    = mymont$SLXCLAY
             ntext      = mymont$NTEXT.SOIL[nzg]

             soil.prop  = soil.params(ntext,isoilflg,slxsand,slxclay)

             #----- Mean diurnal cycle. ---------------------------------------------------#
             dcycmean                = list()
             dcycmean$gpp            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$plresp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$leaf.resp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$root.resp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$growth.resp    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hetresp        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$nep            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$nee            = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$sens           = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxlc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxwc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$hflxgc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$latent         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$et             = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxlc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxwc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wflxgc         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$evap           = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$transp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.temp       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.temp       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$leaf.temp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$wood.temp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$gnd.temp       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.shv        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.shv        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$gnd.shv        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.co2        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.co2        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.prss       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$can.prss       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$atm.vels       = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$fs.open        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rain           = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort.beam    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort.diff    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlong          = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rshort.gnd     = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlong.gnd      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlongup        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$albedo         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$albedo.beam    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$albedo.diff    = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmean$rlong.albedo   = matrix(data=0,nrow=totmon,ncol=ndcycle)

             dcycmsqu             = list()
             dcycmsqu$gpp         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$plresp      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$leaf.resp   = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$root.resp   = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hetresp     = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$nep         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$nee         = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$sens        = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxlc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxwc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$hflxgc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$et          = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$latent      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxlc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxwc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$wflxgc      = matrix(data=0,nrow=totmon,ncol=ndcycle)
             dcycmsqu$transp      = matrix(data=0,nrow=totmon,ncol=ndcycle)

             soil.water           = matrix(data=0,nrow=totmon,ncol=nzg)
             soil.temp            = matrix(data=0,nrow=totmon,ncol=nzg)
             soil.mstpot          = matrix(data=0,nrow=totmon,ncol=nzg)

          }#end if
          #--------------------------------------------------------------------------------#

          #----- Load the total number of patches and cohorts. ----------------------------#
          npat.global = c(npat.global, mymont$NPATCHES.GLOBAL)
          ncoh.global = c(ncoh.global, mymont$NCOHORTS.GLOBAL)
          #--------------------------------------------------------------------------------#



          #----- Load the simple variables. -----------------------------------------------#
          gpp             = c(gpp              ,   mymont$MMEAN.GPP                      )
          plresp          = c(plresp           ,   mymont$MMEAN.PLRESP                   )
          leaf.resp       = c(leaf.resp        ,   mymont$MMEAN.LEAF.RESP                )
          root.resp       = c(root.resp        ,   mymont$MMEAN.ROOT.RESP                )
          growth.resp     = c(growth.resp      ,   mymont$MMEAN.GROWTH.RESP              )
          hetresp         = c(hetresp          ,   mymont$MMEAN.RH                       )
          nep             = c(nep              ,   mymont$MMEAN.NEP                      )
          npp             = c(npp              ,   mymont$MMEAN.GPP - mymont$MMEAN.PLRESP)
          nee             = c(nee              ,   mymont$MMEAN.NEE                      )
          sens            = c(sens             , - mymont$MMEAN.SENSIBLE.AC              )
          hflxlc          = c(hflxlc           ,   mymont$MMEAN.SENSIBLE.LC              )
          hflxwc          = c(hflxwc           ,   mymont$MMEAN.SENSIBLE.WC              )
          hflxgc          = c(hflxgc           ,   mymont$MMEAN.SENSIBLE.GC              )
          latent          = c(latent           , - mymont$MMEAN.VAPOR.AC    * alvl       )
          et              = c(et               , - mymont$MMEAN.VAPOR.AC    * day.sec    )
          wflxlc          = c(wflxlc           ,   mymont$MMEAN.VAPOR.LC    * day.sec    )
          wflxwc          = c(wflxwc           ,   mymont$MMEAN.VAPOR.WC    * day.sec    )
          wflxgc          = c(wflxgc           ,   mymont$MMEAN.VAPOR.GC    * day.sec    )
          evap            = c(evap             ,   mymont$MMEAN.EVAP        * day.sec    )
          transp          = c(transp           ,   mymont$MMEAN.TRANSP      * day.sec    )

          mmsqu.gpp       = c(mmsqu.gpp    ,mymont$MMSQU.GPP                             )
          mmsqu.plresp    = c(mmsqu.plresp ,mymont$MMSQU.PLRESP                          )
          mmsqu.leaf.resp = c(mmsqu.leaf.resp ,mymont$MMSQU.PLRESP                       )
          mmsqu.root.resp = c(mmsqu.root.resp ,mymont$MMSQU.PLRESP                       )
          mmsqu.hetresp   = c(mmsqu.hetresp,mymont$MMSQU.RH                              )
          mmsqu.nee       = c(mmsqu.nee    ,mymont$MMSQU.NEE                             )
          mmsqu.sens      = c(mmsqu.sens   ,mymont$MMSQU.SENSIBLE.AC                     )
          mmsqu.hflxlc    = c(mmsqu.hflxlc ,mymont$MMSQU.SENSIBLE.LC                     )
          mmsqu.hflxwc    = c(mmsqu.hflxwc ,mymont$MMSQU.SENSIBLE.WC                     )
          mmsqu.hflxgc    = c(mmsqu.hflxgc ,mymont$MMSQU.SENSIBLE.GC                     )
          mmsqu.et        = c(mmsqu.et     ,mymont$MMSQU.VAPOR.AC  * day.sec * day.sec   )
          mmsqu.latent    = c(mmsqu.latent ,mymont$MMSQU.VAPOR.AC  * alvl    * alvl      )
          mmsqu.wflxlc    = c(mmsqu.wflxlc ,mymont$MMSQU.VAPOR.LC  * day.sec * day.sec   )
          mmsqu.wflxwc    = c(mmsqu.wflxwc ,mymont$MMSQU.VAPOR.WC  * day.sec * day.sec   )
          mmsqu.wflxgc    = c(mmsqu.wflxgc ,mymont$MMSQU.VAPOR.GC  * day.sec * day.sec   )
          mmsqu.evap      = c(mmsqu.evap   ,mymont$MMSQU.EVAP      * day.sec * day.sec   )
          mmsqu.transp    = c(mmsqu.transp ,mymont$MMSQU.TRANSP    * day.sec * day.sec   )

          atm.vels      = c(atm.vels     ,mymont$MMEAN.ATM.VELS                          )
          atm.prss      = c(atm.prss     ,mymont$MMEAN.ATM.PRSS  * 0.01                  )
          atm.temp      = c(atm.temp     ,mymont$MMEAN.ATM.TEMP  - t00                   )
          atm.shv       = c(atm.shv      ,mymont$MMEAN.ATM.SHV   * kg2g                  )
          atm.co2       = c(atm.co2      ,mymont$MMEAN.ATM.CO2                           )

          can.prss      = c(can.prss     ,mymont$MMEAN.CAN.PRSS  * 0.01                  )
          can.temp      = c(can.temp     ,mymont$MMEAN.CAN.TEMP  - t00                   )
          can.shv       = c(can.shv      ,mymont$MMEAN.CAN.SHV   * kg2g                  )
          can.co2       = c(can.co2      ,mymont$MMEAN.CAN.CO2                           )

          gnd.temp      = c(gnd.temp     ,mymont$MMEAN.GND.TEMP  - t00                   )
          gnd.shv       = c(gnd.shv      ,mymont$MMEAN.GND.SHV   * kg2g                  )

          leaf.temp     = c(leaf.temp    ,mymont$MMEAN.LEAF.TEMP  - t00                  )
          wood.temp     = c(wood.temp    ,mymont$MMEAN.WOOD.TEMP  - t00                  )
          rain          = c(rain         ,mymont$MMEAN.PCPG*ddd  * day.sec               )

          fs.open       = c(fs.open      ,mymont$MMEAN.FS.OPEN                           )
          rshort        = c(rshort       ,mymont$MMEAN.RSHORT                            )
          rshort.beam   = c(rshort.beam  ,mymont$MMEAN.RSHORT - mymont$MMEAN.RSHORT.DIFF )
          rshort.diff   = c(rshort.diff  ,mymont$MMEAN.RSHORT.DIFF                       )
          rlong         = c(rlong        ,mymont$MMEAN.RLONG                             )
          rshort.gnd    = c(rshort.gnd   ,mymont$MMEAN.RSHORT.GND                        )
          rlong.gnd     = c(rlong.gnd    ,mymont$MMEAN.RLONG.GND                         )
          rlongup       = c(rlongup      ,mymont$MMEAN.RLONGUP                           )
          albedo        = c(albedo       ,mymont$MMEAN.ALBEDO                            )
          albedo.beam   = c(albedo.beam  ,mymont$MMEAN.ALBEDO.BEAM                       )
          albedo.diff   = c(albedo.diff  ,mymont$MMEAN.ALBEDO.DIFFUSE                    )
          rlong.albedo  = c(rlong.albedo ,mymont$MMEAN.RLONG.ALBEDO                      )
          #--------------------------------------------------------------------------------#



          #------ Read in soil properties. ------------------------------------------------#
          soil.temp  [m,] =   mymont$MMEAN.SOIL.TEMP - t00
          soil.water [m,] =   mymont$MMEAN.SOIL.WATER
          soil.mstpot[m,] = - mymont$MMEAN.SOIL.MSTPOT
          #--------------------------------------------------------------------------------#



          #----- Read workload, and retrieve only the current month. ----------------------#
          workload   = c(workload  , mymont$WORKLOAD[month])
          specwork   = c(specwork  , mymont$WORKLOAD[month]/sum(mymont$SIPA.N,na.rm=TRUE))
          #--------------------------------------------------------------------------------#



          #----- Disturbance rates. -------------------------------------------------------#
          dist  [m,1:nlu,1:nlu] = mymont$DISTURBANCE.RATES[,1:nlu,1:nlu]
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #       Read the mean diurnal cycle and the mean sum of the squares.             #
          #--------------------------------------------------------------------------------#
          dcycmean$gpp         [m,] = mymont$QMEAN.GPP
          dcycmean$plresp      [m,] = mymont$QMEAN.PLRESP
          dcycmean$leaf.resp   [m,] = mymont$QMEAN.LEAF.RESP
          dcycmean$root.resp   [m,] = mymont$QMEAN.ROOT.RESP
          dcycmean$hetresp     [m,] = mymont$QMEAN.RH
          dcycmean$nep         [m,] = mymont$QMEAN.NEP
          dcycmean$nee         [m,] = mymont$QMEAN.NEE
          dcycmean$sens        [m,] = - mymont$QMEAN.SENSIBLE.AC
          dcycmean$hflxlc      [m,] = mymont$QMEAN.SENSIBLE.LC
          dcycmean$hflxwc      [m,] = mymont$QMEAN.SENSIBLE.WC
          dcycmean$hflxgc      [m,] = mymont$QMEAN.SENSIBLE.GC
          dcycmean$et          [m,] = - mymont$QMEAN.VAPOR.AC         * day.sec
          dcycmean$latent      [m,] = - mymont$QMEAN.VAPOR.AC         * alvl
          dcycmean$wflxlc      [m,] = mymont$QMEAN.VAPOR.LC           * day.sec
          dcycmean$wflxwc      [m,] = mymont$QMEAN.VAPOR.WC           * day.sec
          dcycmean$wflxgc      [m,] = mymont$QMEAN.VAPOR.GC           * day.sec
          dcycmean$evap        [m,] = ( mymont$QMEAN.VAPOR.GC
                                      + mymont$QMEAN.VAPOR.WC
                                      + mymont$QMEAN.VAPOR.LC )       * day.sec
          dcycmean$transp      [m,] = mymont$QMEAN.TRANSP             * day.sec
          dcycmean$atm.temp    [m,] = mymont$QMEAN.ATM.TEMP           - t00
          dcycmean$can.temp    [m,] = mymont$QMEAN.CAN.TEMP           - t00
          dcycmean$leaf.temp   [m,] = mymont$QMEAN.LEAF.TEMP          - t00
          dcycmean$wood.temp   [m,] = mymont$QMEAN.WOOD.TEMP          - t00
          dcycmean$gnd.temp    [m,] = mymont$QMEAN.GND.TEMP           - t00
          dcycmean$atm.shv     [m,] = mymont$QMEAN.ATM.SHV            * kg2g
          dcycmean$can.shv     [m,] = mymont$QMEAN.CAN.SHV            * kg2g
          dcycmean$gnd.shv     [m,] = mymont$QMEAN.GND.SHV            * kg2g
          dcycmean$atm.co2     [m,] = mymont$QMEAN.ATM.CO2
          dcycmean$can.co2     [m,] = mymont$QMEAN.CAN.CO2
          dcycmean$atm.vels    [m,] = mymont$QMEAN.ATM.VELS
          dcycmean$atm.prss    [m,] = mymont$QMEAN.ATM.PRSS * 0.01
          dcycmean$can.prss    [m,] = mymont$QMEAN.CAN.PRSS * 0.01
          dcycmean$fs.open     [m,] = mymont$QMEAN.FS.OPEN
          dcycmean$rain        [m,] = mymont$QMEAN.PCPG               * day.sec
          dcycmean$rshort      [m,] = mymont$QMEAN.RSHORT
          dcycmean$rshort.beam [m,] = mymont$QMEAN.RSHORT - mymont$QMEAN.RSHORT.DIFF
          dcycmean$rshort.diff [m,] = mymont$QMEAN.RSHORT.DIFF
          dcycmean$rlong       [m,] = mymont$QMEAN.RLONG
          dcycmean$rshort.gnd  [m,] = mymont$QMEAN.RSHORT.GND
          dcycmean$rlong.gnd   [m,] = mymont$QMEAN.RLONG.GND
          dcycmean$rlongup     [m,] = mymont$QMEAN.RLONGUP
          dcycmean$albedo      [m,] = mymont$QMEAN.ALBEDO
          dcycmean$albedo.beam [m,] = mymont$QMEAN.ALBEDO.BEAM
          dcycmean$albedo.diff [m,] = mymont$QMEAN.ALBEDO.DIFFUSE
          dcycmean$rlong.albedo[m,] = mymont$QMEAN.RLONG.ALBEDO

          dcycmsqu$gpp         [m,] = mymont$QMSQU.GPP
          dcycmsqu$plresp      [m,] = mymont$QMSQU.PLRESP
          dcycmsqu$leaf.resp   [m,] = mymont$QMSQU.LEAF.RESP
          dcycmsqu$root.resp   [m,] = mymont$QMSQU.ROOT.RESP
          dcycmsqu$hetresp     [m,] = mymont$QMSQU.RH
          dcycmsqu$nep         [m,] = mymont$QMSQU.NEP
          dcycmsqu$nee         [m,] = mymont$QMSQU.NEE
          dcycmsqu$sens        [m,] = mymont$QMSQU.SENSIBLE.AC
          dcycmsqu$hflxlc      [m,] = mymont$QMSQU.SENSIBLE.LC
          dcycmsqu$hflxwc      [m,] = mymont$QMSQU.SENSIBLE.WC
          dcycmsqu$hflxgc      [m,] = mymont$QMSQU.SENSIBLE.GC
          dcycmsqu$et          [m,] = mymont$QMSQU.VAPOR.AC    * day.sec * day.sec
          dcycmsqu$latent      [m,] = mymont$QMSQU.VAPOR.AC    * alvl    * alvl
          dcycmsqu$wflxlc      [m,] = mymont$QMSQU.VAPOR.WC    * day.sec * day.sec
          dcycmsqu$wflxwc      [m,] = mymont$QMSQU.VAPOR.LC    * day.sec * day.sec
          dcycmsqu$wflxgc      [m,] = mymont$QMSQU.VAPOR.GC    * day.sec * day.sec
          dcycmsqu$transp      [m,] = mymont$QMSQU.TRANSP      * day.sec * day.sec
          #--------------------------------------------------------------------------------#


          #---- Read in the site-level area. ----------------------------------------------#
          areasi     = mymont$AREA.SI
          npatches   = mymont$SIPA.N
          #--------------------------------------------------------------------------------#


          #----- Read a few patch-level variables. ----------------------------------------#
          areapa     = mymont$AREA * rep(areasi,times=npatches)
          lupa       = mymont$DIST.TYPE
          agepa      = mymont$AGE
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #    If this is a biomass initialisation, or a run with anthropogenic            #
          # disturbance, we must jitter the age so we can distinguish the patches.         #
          #--------------------------------------------------------------------------------#
          sameage        = duplicated(agepa)
          agepa[sameage] = jitter(x=agepa[sameage],amount=0.4)
          #--------------------------------------------------------------------------------#


          #--------------------------------------------------------------------------------#
          #     Read the cohort-level variables.  Because empty patchs do exist (deserts), #
          # we must check whether there is any cohort to be read.  If not, assign NA to    #
          # all variables.                                                                 #
          #--------------------------------------------------------------------------------#
          ncohorts   = mymont$PACO.N
          if (any (ncohorts > 0)){
             areaconow  = rep(areapa,times=ncohorts)

             #----- Define the land use classes. ------------------------------------------#
             luconow    = rep(lupa,times=ncohorts)

             #----- Define the DBH classes. -----------------------------------------------#
             dbhconow   = mymont$DBH
             dbhbks     = c(-Inf,seq(from=ddbh,to=(ndbh-1)*ddbh,by=ddbh),Inf)
             dbhcut     = cut(dbhconow,breaks=dbhbks)
             dbhlevs    = levels(dbhcut)
             dbhfac     = match(dbhcut,dbhlevs)

             #----- Define the age classes. -----------------------------------------------#
             ageconow   = rep(x=agepa,times=ncohorts)
             agebks     = c(-Inf,seq(from=dage,to=(nage-1)*dage,by=dage),Inf)
             agecut     = cut(ageconow,breaks=agebks)
             agelevs    = levels(agecut)
             agefac     = match(agecut,agelevs)

             agepacut   = cut(agepa,breaks=agebks)
             agepafac   = match(agepacut,agelevs)
             areapaage  = tapply(X=areapa,INDEX=agepafac,sum,na.rm=TRUE)
             areaage    = areapaage[as.character(agefac)]

             pftconow        = mymont$PFT
             nplantconow     = mymont$NPLANT
             heightconow     = mymont$HITE
             baconow         = mymont$BA.CO
             agbconow        = mymont$AGB.CO
             bseedsconow     = mymont$BSEEDS.CO
             laiconow        = mymont$LAI.CO
             waiconow        = mymont$WAI.CO
             taiconow        = laiconow + waiconow
             gppconow        = mymont$MMEAN.GPP.CO
             leafrespconow   = mymont$MMEAN.LEAF.RESP.CO
             rootrespconow   = mymont$MMEAN.ROOT.RESP.CO
             growthrespconow = mymont$MMEAN.GROWTH.RESP.CO
             respconow       = ( mymont$MMEAN.LEAF.RESP.CO   + mymont$MMEAN.ROOT.RESP.CO
                               + mymont$MMEAN.GROWTH.RESP.CO + mymont$MMEAN.STORAGE.RESP.CO
                               + mymont$MMEAN.VLEAF.RESP.CO  )
             nppconow        = gppconow-respconow
             cbalconow       = mymont$MMEAN.CB
             mcostconow      = ( mymont$MMEAN.LEAF.MAINTENANCE
                               + mymont$MMEAN.ROOT.MAINTENANCE )
             ldropconow      = mymont$MMEAN.LEAF.DROP.CO
             cbrbarconow     = mymont$CBR.BAR
             ncbmortconow    = mymont$MMEAN.MORT.RATE[,2]
             fsoconow        = mymont$MMEAN.FS.OPEN.CO
             lightconow      = mymont$MMEAN.LIGHT.LEVEL
             lambdaconow     = mymont$MMEAN.LAMBDA.LIGHT.CO
             beamextconow    = mymont$MMEAN.BEAMEXT.LEVEL
             diffextconow    = mymont$MMEAN.BEAMEXT.LEVEL
             parlconow       = mymont$MMEAN.PAR.L

             baliveconow     = mymont$BALIVE
             bdeadconow      = mymont$BDEAD
             bleafconow      = mymont$BLEAF
             brootconow      = mymont$BROOT
             bswoodconow     = mymont$BSAPWOOD
             bstoreconow     = mymont$BSTORAGE


             sel               = laiconow > 1.e-10
             demandconow       = mymont$MMEAN.PSI.OPEN     * day.sec
             supplyconow       = mymont$MMEAN.WATER.SUPPLY * day.sec
             gpplconow         = gppconow
             demandconow[sel]  = demandconow[sel] / laiconow[sel]
             supplyconow[sel]  = supplyconow[sel] / laiconow[sel]
             gpplconow  [sel]  = nplantconow[sel] * gppconow[sel] / laiconow[sel]
             demandconow[!sel] = 0.
             supplyconow[!sel] = 0.
             gpplconow  [!sel] = 0.
             #-----------------------------------------------------------------------------#
          }else{
             #----- Make everything NA. ---------------------------------------------------#
             areaconow       = NA
             luconow         = NA
             dbhconow        = NA
             dbhbks          = NA
             dbhcut          = NA
             dbhlevs         = NA
             dbhfac          = NA
             ageconow        = NA
             agebks          = NA
             agecut          = NA
             agelevs         = NA
             agefac          = NA
             agepacut        = NA
             agepafac        = NA
             areapaage       = NA
             areaage         = NA
             pftconow        = NA
             nplantconow     = NA
             heightconow     = NA
             agbconow        = NA
             baconow         = NA
             bseedsconow     = NA
             laiconow        = NA
             waiconow        = NA
             taiconow        = NA
             gppconow        = NA
             gpplconow       = NA
             leafrespconow   = NA
             rootrespconow   = NA
             growthrespconow = NA
             respconow       = NA 
             nppconow        = NA 
             cbalconow       = NA 
             mcostconow      = NA 
             ldropconow      = NA 
             cbrbarconow     = NA 
             ncbmortconow    = NA 
             fsoconow        = NA 
             lightconow      = NA 
             lambdaconow     = NA 
             beamextconow    = NA 
             diffextconow    = NA 
             parlconow       = NA 
             demandconow     = NA 
             supplyconow     = NA 
             baliveconow     = NA 
             bdeadconow      = NA 
             bleafconow      = NA 
             brootconow      = NA 
             bswoodconow     = NA 
             bstoreconow     = NA 
          }#end if

          #----- Define some classes that can be defined even when no cohorts exist. ------#
          agebks     = c(-Inf,seq(from=dage,to=(nage-1)*dage,by=dage),Inf)
          agepacut   = cut(agepa,breaks=agebks)
          agepafac   = match(agepacut,agelevs)
          areapaage  = tapply(X=areapa,INDEX=agepafac,sum,na.rm=TRUE)
          areaage    = areapaage[as.character(agefac)]



          #--------------------------------------------------------------------------------#
          #     Build the PFT arrays.                                                      #
          #--------------------------------------------------------------------------------#
          for (p in 1:npft){
              if (all(is.na(pftconow))){
                 sel      = rep(FALSE,times=length(pftconow))
              }else{
                 sel      = pftconow == p
              }#end if
              if (any(sel)){
                 laipft    [m,p] = laipft     [m,p] + sum(laiconow [sel] * areaconow[sel])
                 waipft    [m,p] = waipft     [m,p] + sum(waiconow [sel] * areaconow[sel])
                 taipft    [m,p] = taipft     [m,p] + sum(taiconow [sel] * areaconow[sel])

                 basareapft[m,p]    = ( basareapft [m,p] 
                                      + sum( nplantconow[sel] * baconow [sel]   
                                           * areaconow[sel]))
                 agbpft    [m,p]    = ( agbpft [m,p] 
                                      + sum( nplantconow[sel] * agbconow [sel]   
                                           * areaconow[sel]))
                 bseedspft [m,p]    = ( bseedspft [m,p]
                                      + sum( nplantconow[sel] * bseedsconow [sel] 
                                           * areaconow [sel]))
                 gpppft    [m,p]    = ( gpppft [m,p]
                                      + sum( nplantconow[sel] * gppconow [sel]
                                           * areaconow[sel]))
                 npppft    [m,p]    = ( npppft [m,p]
                                      + sum( nplantconow[sel] * nppconow [sel]  
                                           * areaconow[sel]))
                 mcopft    [m,p]    = ( mcopft [m,p]
                                      + sum( nplantconow[sel] * mcostconow [sel]
                                           * areaconow[sel]))
                 cbapft    [m,p]    = ( cbapft [m,p] 
                                      + sum( nplantconow[sel] * cbalconow [sel]  
                                           * areaconow[sel]))
                 ldroppft  [m,p]    = ( ldroppft [m,p] 
                                      + sum( nplantconow[sel] * ldropconow [sel]  
                                           * areaconow[sel]))
                 balivepft [m,p]    = ( balivepft [m,p]
                                      + sum( nplantconow[sel] * baliveconow[sel]
                                           * areaconow[sel]))
                 bdeadpft  [m,p]    = ( bdeadpft [m,p]
                                      + sum( nplantconow[sel] * bdeadconow[sel]
                                           * areaconow[sel]))
                 bleafpft  [m,p]    = ( bleafpft [m,p]
                                      + sum( nplantconow[sel] * bleafconow[sel]
                                           * areaconow[sel]))
                 brootpft  [m,p]    = ( brootpft [m,p]
                                      + sum( nplantconow[sel] * brootconow[sel]
                                           * areaconow[sel]))
                 bswoodpft [m,p]    = ( bswoodpft [m,p]
                                      + sum( nplantconow[sel] * bswoodconow[sel]
                                           * areaconow[sel]))
                 bstorepft [m,p]    = ( bstorepft [m,p]
                                      + sum( nplantconow[sel] * bstoreconow[sel]
                                           * areaconow[sel]))
                 leafresppft[m,p]   = ( leafresppft [m,p] 
                                      + sum( nplantconow[sel] * leafrespconow[sel]
                                           * areaconow[sel]))
                 rootresppft[m,p]   = ( rootresppft [m,p] 
                                      + sum( nplantconow[sel] * rootrespconow[sel]
                                           * areaconow[sel]))
                 growthresppft[m,p] = ( growthresppft [m,p] 
                                      + sum( nplantconow[sel] * growthrespconow[sel]
                                           * areaconow[sel]))
              }
          }
          #------ Find the total. ---------------------------------------------------------#
          laipft       [m,npft+1] = sum(laipft       [m,1:npft],na.rm=TRUE)
          waipft       [m,npft+1] = sum(waipft       [m,1:npft],na.rm=TRUE)
          taipft       [m,npft+1] = sum(taipft       [m,1:npft],na.rm=TRUE)
          agbpft       [m,npft+1] = sum(agbpft       [m,1:npft],na.rm=TRUE)
          bseedspft    [m,npft+1] = sum(bseedspft    [m,1:npft],na.rm=TRUE)
          gpppft       [m,npft+1] = sum(gpppft       [m,1:npft],na.rm=TRUE)
          npppft       [m,npft+1] = sum(npppft       [m,1:npft],na.rm=TRUE)
          mcopft       [m,npft+1] = sum(mcopft       [m,1:npft],na.rm=TRUE)
          cbapft       [m,npft+1] = sum(cbapft       [m,1:npft],na.rm=TRUE)
          ldroppft     [m,npft+1] = sum(ldroppft     [m,1:npft],na.rm=TRUE)
          balivepft    [m,npft+1] = sum(balivepft    [m,1:npft],na.rm=TRUE)
          bdeadpft     [m,npft+1] = sum(bdeadpft     [m,1:npft],na.rm=TRUE)
          bleafpft     [m,npft+1] = sum(bleafpft     [m,1:npft],na.rm=TRUE)
          brootpft     [m,npft+1] = sum(brootpft     [m,1:npft],na.rm=TRUE)
          bswoodpft    [m,npft+1] = sum(bswoodpft    [m,1:npft],na.rm=TRUE)
          bstorepft    [m,npft+1] = sum(bstorepft    [m,1:npft],na.rm=TRUE)
          basareapft   [m,npft+1] = sum(basareapft   [m,1:npft],na.rm=TRUE)
          leafresppft  [m,npft+1] = sum(leafresppft  [m,1:npft],na.rm=TRUE)
          rootresppft  [m,npft+1] = sum(rootresppft  [m,1:npft],na.rm=TRUE)
          growthresppft[m,npft+1] = sum(growthresppft[m,1:npft],na.rm=TRUE)
          #--------------------------------------------------------------------------------#




          #--------------------------------------------------------------------------------#
          #     Build the LU arrays.                                                       #
          #--------------------------------------------------------------------------------#
          for (l in 1:nlu){
             selpa    = lupa    == l
             if (all(is.na(luconow))){
                sel      = rep(FALSE,times=length(luconow))
             }else{
                sel      = luconow == l
             }#end if
             if (any(sel)){
                lailu    [m,l] = lailu [m,l]    + sum(laiconow [sel] * areaconow[sel])
                basarealu[m,l] = basarealu [m,l] + 
                                 sum(nplantconow[sel] * baconow [sel]    * areaconow[sel])
                agblu [m,l]    = agblu [m,l] + 
                                 sum(nplantconow[sel] * agbconow [sel]   * areaconow[sel])
                gpplu [m,l]    = gpplu [m,l] + 
                                 sum(nplantconow[sel] * gppconow [sel]   * areaconow[sel])
                npplu [m,l]    = npplu [m,l] +
                                 sum(nplantconow[sel] * nppconow [sel]   * areaconow[sel])
             }#end if
             arealu [m,l]   = arealu [m,l] + sum(areapa[selpa])
          }#end for
          #------ Find the total. ---------------------------------------------------------#
          lailu    [m,nlu+1] = sum(lailu    [m,1:nlu], na.rm=TRUE)
          basarealu[m,nlu+1] = sum(basarealu[m,1:nlu], na.rm=TRUE)
          agblu    [m,nlu+1] = sum(agblu    [m,1:nlu], na.rm=TRUE)
          gpplu    [m,nlu+1] = sum(gpplu    [m,1:nlu], na.rm=TRUE)
          npplu    [m,nlu+1] = sum(npplu    [m,1:nlu], na.rm=TRUE)
          arealu   [m,nlu+1] = sum(arealu   [m,1:nlu], na.rm=TRUE)
          #--------------------------------------------------------------------------------#




          #--------------------------------------------------------------------------------#
          #     Build the size (DBH) structure arrays.                                     #
          #--------------------------------------------------------------------------------#
          for (d in 1:ndbh){
             if (all(is.na(dbhfac))){
                seldbh  = rep(FALSE,times=length(dbhfac))
             }else{
                seldbh  = dbhfac == d
             }#end if
             for (p in 1:npft){
                 selpft   = pftconow == p
                 sel      = selpft & seldbh
                 if (any(sel)){
                    laipftdbh [m,d,p] = laipftdbh [m,d,p] + 
                                        sum(laiconow [sel] * areaconow[sel])
                    waipftdbh [m,d,p] = waipftdbh [m,d,p] + 
                                        sum(waiconow [sel] * areaconow[sel])
                    taipftdbh [m,d,p] = taipftdbh [m,d,p] + 
                                        sum(taiconow [sel] * areaconow[sel])
                    agbpftdbh [m,d,p] = agbpftdbh [m,d,p] + 
                                        sum( nplantconow[sel] * agbconow   [sel]
                                           * areaconow[sel])
                    gpppftdbh [m,d,p] = gpppftdbh [m,d,p] + 
                                        sum( nplantconow[sel] * gppconow   [sel]
                                           * areaconow[sel])
                    npppftdbh [m,d,p] = npppftdbh [m,d,p] +
                                        sum( nplantconow[sel] * nppconow   [sel]
                                           * areaconow[sel])
                    mcopftdbh [m,d,p] = mcopftdbh [m,d,p] +
                                        sum( nplantconow[sel] * mcostconow [sel]
                                           * areaconow[sel])
                    cbapftdbh [m,d,p] = cbapftdbh [m,d,p] +
                                        sum( nplantconow[sel] * cbalconow  [sel]
                                           * areaconow[sel])
                 }
             }
          }
          #--------------------------------------------------------------------------------#




          #--------------------------------------------------------------------------------#
          #      Build the age structure arrays.                                           #
          #--------------------------------------------------------------------------------#
          for (a in 1:nage){
             if (all(is.na(agefac))){
                selage   = rep(FALSE,times=length(agefac))
             }else{
                selage  = agefac == a
             }#end if
             for (p in 1:npft){
                if (all(is.na(pftconow))){
                   selpft   = rep(FALSE,times=length(pftconow))
                }else{
                   selpft   = pftconow == p
                }#end if
                sel    = selpft & selage
                if (any(sel)){
                    laipftage [m,a,p] = laipftage [m,a,p]  + 
                                        sum(laiconow [sel] * areaconow[sel] / areaage[sel])
                    waipftage [m,a,p] = waipftage [m,a,p]  + 
                                        sum(waiconow [sel] * areaconow[sel] / areaage[sel])
                    taipftage [m,a,p] = taipftage [m,a,p]  + 
                                        sum(taiconow [sel] * areaconow[sel] / areaage[sel])
                    agbpftage [m,a,p] = agbpftage [m,a,p]  + 
                                        sum( nplantconow[sel] * agbconow [sel]   
                                           * areaconow[sel]   / areaage[sel])
                    gpppftage [m,a,p] = gpppftage [m,a,p]  + 
                                        sum( nplantconow[sel] * gppconow [sel]   
                                           * areaconow[sel]   / areaage[sel])
                    npppftage [m,a,p] = npppftage [m,a,p]  + 
                                        sum( nplantconow[sel] * nppconow [sel]
                                           * areaconow[sel]   / areaage[sel])
                    mcopftage [m,a,p] = mcopftage [m,a,p]  + 
                                        sum( nplantconow[sel] * mcostconow [sel]
                                           * areaconow[sel]   / areaage[sel])
                    cbapftage [m,a,p] = cbapftage [m,a,p]  + 
                                        sum( nplantconow[sel] * cbalconow [sel]
                                           * areaconow[sel]   / areaage[sel])
                }
             }
          }
          #--------------------------------------------------------------------------------#




          #--------------------------------------------------------------------------------#
          #       Build the derived variables.                                             #
          #--------------------------------------------------------------------------------#
          npp   = c(npp  ,sum(npppft  [m,]) )
          mco   = c(mco  ,sum(mcopft  [m,]) )
          cba   = c(cba  ,sum(cbapft  [m,]) )
          lai   = c(lai  ,sum(laipft  [m,]) )
          wai   = c(wai  ,sum(waipft  [m,]) )
          tai   = c(tai  ,sum(taipft  [m,]) )
          agb   = c(agb  ,sum(agbpft  [m,]) )
          ldrop = c(ldrop,sum(ldroppft[m,]) )
          #--------------------------------------------------------------------------------#




          #--------------------------------------------------------------------------------#
          #      Build the cohort-level lists if this is the right month.                  #
          #--------------------------------------------------------------------------------#
          if (month %in% sasmonth){
             cyear  = substring(10000 + year,2,5)
             cmonth = substring(100+month,2,3)
             labwhen     = paste("y",cyear,"m",cmonth,sep="")
             #----- Binding the current cohorts. ------------------------------------------#
             lightco  [[labwhen]] = lightconow
             beamextco[[labwhen]] = beamextconow
             diffextco[[labwhen]] = diffextconow
             parlco   [[labwhen]] = parlconow
             lambdaco [[labwhen]] = lambdaconow
             gppco    [[labwhen]] = gppconow
             gpplco   [[labwhen]] = gpplconow
             respco   [[labwhen]] = respconow
             nppco    [[labwhen]] = nppconow
             cbrbarco [[labwhen]] = cbrbarconow
             cbalco   [[labwhen]] = cbalconow
             mcostco  [[labwhen]] = mcostconow
             ncbmortco[[labwhen]] = ncbmortconow
             agbco    [[labwhen]] = agbconow
             fsoco    [[labwhen]] = fsoconow
             nplantco [[labwhen]] = nplantconow * areaconow
             heightco [[labwhen]] = heightconow
             baco     [[labwhen]] = nplantconow * baconow * areaconow
             pftco    [[labwhen]] = pftconow
             dbhco    [[labwhen]] = dbhconow
             laico    [[labwhen]] = laiconow
             waico    [[labwhen]] = waiconow
             taico    [[labwhen]] = taiconow
             ageco    [[labwhen]] = ageconow
             areaco   [[labwhen]] = areaconow
             demandco [[labwhen]] = demandconow
             supplyco [[labwhen]] = supplyconow
             baliveco [[labwhen]] = baliveconow
             bdeadco  [[labwhen]] = bdeadconow
             bleafco  [[labwhen]] = bleafconow
             brootco  [[labwhen]] = brootconow
             bswoodco [[labwhen]] = bswoodconow
             bstoreco [[labwhen]] = bstoreconow
          } #end if month=sasmonth
          #--------------------------------------------------------------------------------#

      }# end for, month
   }#end for, year



   #----- Find the date variables. --------------------------------------------------------#
   thismonth = chron(thismonth,out.format=c(dates="day-mon-yr",times=NULL))
   mmonth    = months(thismonth)
   mfac      = nummonths(thismonth)
   moncensus = census(pop=mfac,categ=seq(from=1,to=12,by=1))
   moncnt    = matrix(data=rep(moncensus,times=ndcycle),nrow=12)


   #---------------------------------------------------------------------------------------#
   #      Here we find the monthly means for month, then compute the standard deviation.   #
   #---------------------------------------------------------------------------------------#
   print ("    - Finding the monthly mean...")
   print ("      * Aggregating the monthly mean...")
   mont12mn             = list()
   mont12mn$gpp         = tapply(X=gpp          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$plresp      = tapply(X=plresp       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$leaf.resp   = tapply(X=leaf.resp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$root.resp   = tapply(X=root.resp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$growth.resp = tapply(X=growth.resp  ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$hetresp     = tapply(X=hetresp      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$nee         = tapply(X=nee          ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$sens        = tapply(X=sens         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$hflxlc      = tapply(X=hflxlc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$hflxwc      = tapply(X=hflxwc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$hflxgc      = tapply(X=hflxgc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$et          = tapply(X=et           ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$latent      = tapply(X=latent       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$wflxlc      = tapply(X=wflxlc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$wflxwc      = tapply(X=wflxwc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$wflxgc      = tapply(X=wflxgc       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$evap        = tapply(X=evap         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$transp      = tapply(X=transp       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$rain        = tapply(X=rain         ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$atm.temp    = tapply(X=atm.temp     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$rshort      = tapply(X=rshort       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$rlong       = tapply(X=rlong        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$atm.shv     = tapply(X=atm.shv      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$atm.co2     = tapply(X=atm.co2      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$atm.prss    = tapply(X=atm.prss     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$atm.vels    = tapply(X=atm.vels     ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12mn$soil.temp   = qapply(mat=soil.temp  ,index=mfac,bycol=T,func=mean,na.rm=T)
   mont12mn$soil.water  = qapply(mat=soil.water ,index=mfac,bycol=T,func=mean,na.rm=T)
   mont12mn$soil.mstpot = qapply(mat=soil.mstpot,index=mfac,bycol=T,func=mean,na.rm=T)
   #----- Find the mean sum of squares. ---------------------------------------------------#
   print ("      * Aggregating the monthly mean sum of squares...")
   mont12sq           = list()
   mont12sq$gpp       = tapply(X=mmsqu.gpp       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$plresp    = tapply(X=mmsqu.plresp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$leaf.resp = tapply(X=mmsqu.leaf.resp ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$root.resp = tapply(X=mmsqu.root.resp ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$hetresp   = tapply(X=mmsqu.hetresp   ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$nee       = tapply(X=mmsqu.nee       ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$sens      = tapply(X=mmsqu.sens      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$hflxlc    = tapply(X=mmsqu.hflxlc    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$hflxwc    = tapply(X=mmsqu.hflxwc    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$hflxgc    = tapply(X=mmsqu.hflxgc    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$et        = tapply(X=mmsqu.et        ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$latent    = tapply(X=mmsqu.latent    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$wflxlc    = tapply(X=mmsqu.wflxlc    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$wflxwc    = tapply(X=mmsqu.wflxwc    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$wflxgc    = tapply(X=mmsqu.wflxgc    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$evap      = tapply(X=mmsqu.evap      ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   mont12sq$transp    = tapply(X=mmsqu.transp    ,INDEX=mfac,FUN=mean,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#
   #   Here we convert the sum of squares into standard deviation. The standard devi-      #
   # ation can be written in two different ways, and we will use the latter because it     #
   # doesn't require previous knowledge of the mean.                                       #
   #              __________________          ____________________________________         #
   #             / SUM_i[X_i - Xm]²          /  / SUM_i[X_i²]        \      1              #
   # sigma = \  /  ----------------   =  \  /  |  -----------  - Xm²  | ---------          #
   #          \/       N - 1              \/    \      N             /   1 - 1/N           #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   print ("      * Finding the standard deviation...")
   srnorm1 = sqrt(1./(1. - 1. / moncensus))
   srnorm1[!is.finite(srnorm1)] = 0.
   mont12sd            = list()
   mont12sd$gpp        = sqrt(mont12sq$gpp        - mont12mn$gpp^2        ) * srnorm1
   mont12sd$plresp     = sqrt(mont12sq$plresp     - mont12mn$plresp^2     ) * srnorm1
   mont12sd$leaf.resp  = sqrt(mont12sq$leaf.resp  - mont12mn$leaf.resp^2  ) * srnorm1
   mont12sd$root.resp  = sqrt(mont12sq$root.resp  - mont12mn$root.resp^2  ) * srnorm1
   mont12sd$hetresp    = sqrt(mont12sq$hetresp    - mont12mn$hetresp^2    ) * srnorm1
   mont12sd$nee        = sqrt(mont12sq$nee        - mont12mn$nee^2        ) * srnorm1
   mont12sd$sens       = sqrt(mont12sq$sens       - mont12mn$sens^2       ) * srnorm1
   mont12sd$hflxlc     = sqrt(mont12sq$hflxlc     - mont12mn$hflxlc^2     ) * srnorm1
   mont12sd$hflxwc     = sqrt(mont12sq$hflxwc     - mont12mn$hflxwc^2     ) * srnorm1
   mont12sd$hflxgc     = sqrt(mont12sq$hflxgc     - mont12mn$hflxgc^2     ) * srnorm1
   mont12sd$et         = sqrt(mont12sq$et         - mont12mn$et^2         ) * srnorm1
   mont12sd$latent     = sqrt(mont12sq$latent     - mont12mn$latent^2     ) * srnorm1
   mont12sd$wflxlc     = sqrt(mont12sq$wflxlc     - mont12mn$wflxlc^2     ) * srnorm1
   mont12sd$wflxwc     = sqrt(mont12sq$wflxwc     - mont12mn$wflxwc^2     ) * srnorm1
   mont12sd$wflxgc     = sqrt(mont12sq$wflxgc     - mont12mn$wflxgc^2     ) * srnorm1
   mont12sd$evap       = sqrt(mont12sq$evap       - mont12mn$evap^2       ) * srnorm1
   mont12sd$transp     = sqrt(mont12sq$transp     - mont12mn$transp^2     ) * srnorm1
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   mont12sd$gpp        [!is.finite(mont12mn$gpp        )] = 0.
   mont12sd$plresp     [!is.finite(mont12mn$plresp     )] = 0.
   mont12sd$leaf.resp  [!is.finite(mont12mn$leaf.resp  )] = 0.
   mont12sd$root.resp  [!is.finite(mont12mn$root.resp  )] = 0.
   mont12sd$hetresp    [!is.finite(mont12mn$hetresp    )] = 0.
   mont12sd$nee        [!is.finite(mont12mn$nee        )] = 0.
   mont12sd$sens       [!is.finite(mont12mn$sens       )] = 0.
   mont12sd$hflxlc     [!is.finite(mont12mn$hflxlc     )] = 0.
   mont12sd$hflxlc     [!is.finite(mont12mn$hflxwc     )] = 0.
   mont12sd$hflxgc     [!is.finite(mont12mn$hflxgc     )] = 0.
   mont12sd$et         [!is.finite(mont12mn$et         )] = 0.
   mont12sd$latent     [!is.finite(mont12mn$latent     )] = 0.
   mont12sd$wflxlc     [!is.finite(mont12mn$wflxlc     )] = 0.
   mont12sd$wflxwc     [!is.finite(mont12mn$wflxwc     )] = 0.
   mont12sd$wflxgc     [!is.finite(mont12mn$wflxgc     )] = 0.
   mont12sd$evap       [!is.finite(mont12mn$evap       )] = 0.
   mont12sd$transp     [!is.finite(mont12mn$transp     )] = 0.
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Here we find the Mean diurnal cycle for each month, then compute the standard    #
   # deviation.                                                                            #
   #---------------------------------------------------------------------------------------#
   print ("    - Aggregating the monthly mean of the diurnal cycle...")
   dcyc12mn             = list()
   dcyc12mn$gpp         = qapply(dcycmean$gpp         ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$plresp      = qapply(dcycmean$plresp      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$leaf.resp   = qapply(dcycmean$leaf.resp   ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$root.resp   = qapply(dcycmean$root.resp   ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$hetresp     = qapply(dcycmean$hetresp     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$nep         = qapply(dcycmean$nep         ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$nee         = qapply(dcycmean$nee         ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$sens        = qapply(dcycmean$sens        ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$hflxlc      = qapply(dcycmean$hflxlc      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$hflxwc      = qapply(dcycmean$hflxwc      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$hflxgc      = qapply(dcycmean$hflxgc      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$et          = qapply(dcycmean$et          ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$latent      = qapply(dcycmean$latent      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$wflxlc      = qapply(dcycmean$wflxlc      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$wflxwc      = qapply(dcycmean$wflxwc      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$wflxgc      = qapply(dcycmean$wflxgc      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$evap        = qapply(dcycmean$evap        ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$transp      = qapply(dcycmean$transp      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$atm.temp    = qapply(dcycmean$atm.temp    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$can.temp    = qapply(dcycmean$can.temp    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$leaf.temp   = qapply(dcycmean$leaf.temp   ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$wood.temp   = qapply(dcycmean$wood.temp   ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$gnd.temp    = qapply(dcycmean$gnd.temp    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$atm.shv     = qapply(dcycmean$atm.shv     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$can.shv     = qapply(dcycmean$can.shv     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$gnd.shv     = qapply(dcycmean$gnd.shv     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$atm.co2     = qapply(dcycmean$atm.co2     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$can.co2     = qapply(dcycmean$can.co2     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$atm.prss    = qapply(dcycmean$atm.prss    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$can.prss    = qapply(dcycmean$can.prss    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$atm.vels    = qapply(dcycmean$atm.vels    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$fs.open     = qapply(dcycmean$fs.open     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rain        = qapply(dcycmean$rain        ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rshort      = qapply(dcycmean$rshort      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rshort.beam = qapply(dcycmean$rshort.beam ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rshort.diff = qapply(dcycmean$rshort.diff ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rlong       = qapply(dcycmean$rlong       ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rshort.gnd  = qapply(dcycmean$rshort.gnd  ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rlong.gnd   = qapply(dcycmean$rlong.gnd   ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rlongup     = qapply(dcycmean$rlongup     ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$albedo      = qapply(dcycmean$albedo      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$albedo.beam = qapply(dcycmean$albedo.beam ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$albedo.diff = qapply(dcycmean$albedo.diff ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12mn$rlong.albedo= qapply(dcycmean$rlong.albedo,index=mfac,bycol=T,func=mean,na.rm=T)

   #----- Find the mean sum of squares. ---------------------------------------------------#
   print ("    - Aggregating the monthly mean sum of squares...")
   dcyc12sq            = list()
   dcyc12sq$gpp        = qapply(dcycmsqu$gpp       ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$plresp     = qapply(dcycmsqu$plresp    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$leaf.resp  = qapply(dcycmsqu$leaf.resp ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$root.resp  = qapply(dcycmsqu$root.resp ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$hetresp    = qapply(dcycmsqu$hetresp   ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$nep        = qapply(dcycmsqu$nep       ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$nee        = qapply(dcycmsqu$nee       ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$sens       = qapply(dcycmsqu$sens      ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$hflxlc     = qapply(dcycmsqu$hflxlc    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$hflxwc     = qapply(dcycmsqu$hflxwc    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$hflxgc     = qapply(dcycmsqu$hflxgc    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$et         = qapply(dcycmsqu$et        ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$latent     = qapply(dcycmsqu$latent    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$wflxlc     = qapply(dcycmsqu$wflxlc    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$wflxwc     = qapply(dcycmsqu$wflxwc    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$wflxgc     = qapply(dcycmsqu$wflxgc    ,index=mfac,bycol=T,func=mean,na.rm=T)
   dcyc12sq$transp     = qapply(dcycmsqu$transp    ,index=mfac,bycol=T,func=mean,na.rm=T)

   #---------------------------------------------------------------------------------------#
   #   Here we convert the sum of squares into standard deviation. The standard devi-      #
   # ation can be written in two different ways, and we will use the latter because it     #
   # doesn't require previous knowledge of the mean.                                       #
   #              __________________          ____________________________________         #
   #             / SUM_i[X_i - Xm]²          /  / SUM_i[X_i²]        \      1              #
   # sigma = \  /  ----------------   =  \  /  |  -----------  - Xm²  | ---------          #
   #          \/       N - 1              \/    \      N             /   1 - 1/N           #
   #                                                                                       #
   # srnonm1 is the square root of 1 / (1 - 1/N)                                           #
   #     Find the standard deviation.                                                      #
   #---------------------------------------------------------------------------------------#
   print ("    - Finding the standard deviation...")
   srnorm1 = sqrt(1./(1. - 1. / moncnt))
   srnorm1[!is.finite(srnorm1)] = 0.
   dcyc12sd            = list()
   dcyc12sd$gpp        = sqrt(dcyc12sq$gpp      -dcyc12mn$gpp^2             )*srnorm1
   dcyc12sd$plresp     = sqrt(dcyc12sq$plresp   -dcyc12mn$plresp^2          )*srnorm1
   dcyc12sd$leaf.resp  = sqrt(dcyc12sq$leaf.resp-dcyc12mn$leaf.resp^2       )*srnorm1
   dcyc12sd$root.resp  = sqrt(dcyc12sq$root.resp-dcyc12mn$root.resp^2       )*srnorm1
   dcyc12sd$hetresp    = sqrt(dcyc12sq$hetresp  -dcyc12mn$hetresp^2         )*srnorm1
   dcyc12sd$nep        = sqrt(dcyc12sq$nep      -dcyc12mn$nep^2             )*srnorm1
   dcyc12sd$nee        = sqrt(dcyc12sq$nee      -dcyc12mn$nee^2             )*srnorm1
   dcyc12sd$sens       = sqrt(dcyc12sq$sens     -dcyc12mn$sens^2            )*srnorm1
   dcyc12sd$hflxlc     = sqrt(dcyc12sq$hflxlc   -dcyc12mn$hflxlc^2          )*srnorm1
   dcyc12sd$hflxwc     = sqrt(dcyc12sq$hflxwc   -dcyc12mn$hflxwc^2          )*srnorm1
   dcyc12sd$hflxgc     = sqrt(dcyc12sq$hflxgc   -dcyc12mn$hflxgc^2          )*srnorm1
   dcyc12sd$et         = sqrt(dcyc12sq$et       -dcyc12mn$et^2              )*srnorm1
   dcyc12sd$latent     = sqrt(dcyc12sq$latent   -dcyc12mn$latent^2          )*srnorm1
   dcyc12sd$wflxlc     = sqrt(dcyc12sq$wflxlc   -dcyc12mn$wflxlc^2          )*srnorm1
   dcyc12sd$wflxwc     = sqrt(dcyc12sq$wflxwc   -dcyc12mn$wflxwc^2          )*srnorm1
   dcyc12sd$wflxgc     = sqrt(dcyc12sq$wflxgc   -dcyc12mn$wflxgc^2          )*srnorm1
   dcyc12sd$transp     = sqrt(dcyc12sq$transp   -dcyc12mn$transp^2          )*srnorm1
   #---------------------------------------------------------------------------------------#
   #     Set standard deviations that became NaN to 0.  This usually happens when we run   #
   # the post-processing script when the simulation hasn't run for more than 2 years.  We  #
   # can't find the standard deviation because the number of degrees of freedom becomes 0. #
   #---------------------------------------------------------------------------------------#
   dcyc12sd$gpp        [!is.finite(dcyc12sd$gpp       )] = 0.
   dcyc12sd$plresp     [!is.finite(dcyc12sd$plresp    )] = 0.
   dcyc12sd$leaf.resp  [!is.finite(dcyc12sd$leaf.resp )] = 0.
   dcyc12sd$root.resp  [!is.finite(dcyc12sd$root.resp )] = 0.
   dcyc12sd$hetresp    [!is.finite(dcyc12sd$hetresp   )] = 0.
   dcyc12sd$nep        [!is.finite(dcyc12sd$nep       )] = 0.
   dcyc12sd$nee        [!is.finite(dcyc12sd$nee       )] = 0.
   dcyc12sd$sens       [!is.finite(dcyc12sd$sens      )] = 0.
   dcyc12sd$hflxlc     [!is.finite(dcyc12sd$hflxlc    )] = 0.
   dcyc12sd$hflxwc     [!is.finite(dcyc12sd$hflxwc    )] = 0.
   dcyc12sd$hflxgc     [!is.finite(dcyc12sd$hflxgc    )] = 0.
   dcyc12sd$et         [!is.finite(dcyc12sd$et        )] = 0.
   dcyc12sd$latent     [!is.finite(dcyc12sd$latent    )] = 0.
   dcyc12sd$wflxlc     [!is.finite(dcyc12sd$wflxlc    )] = 0.
   dcyc12sd$wflxwc     [!is.finite(dcyc12sd$wflxwc    )] = 0.
   dcyc12sd$wflxgc     [!is.finite(dcyc12sd$wflxgc    )] = 0.
   dcyc12sd$transp     [!is.finite(dcyc12sd$transp    )] = 0.
   #---------------------------------------------------------------------------------------#



   #----- Find which PFTs, land uses and transitions we need to consider ------------------#
   pftave  = colMeans(agbpft,na.rm=TRUE)
   luave   = colMeans(agblu ,na.rm=TRUE)
   distave = matrix(NA,nrow=3,ncol=3)
   for (jlu in 1:nlu){
      for (ilu in 1:nlu){
          distave[ilu,jlu] = mean(dist[,ilu,jlu])
      }#end for
   }#end for
   selpft  = pftave  > 0.
   sellu   = luave   > 0.
   seldist = distave > 0.


   #----- Determine the last data available. ----------------------------------------------#
   tlast = length(thismonth)

   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for those time series that uses thismonth...             #
   #---------------------------------------------------------------------------------------#
   whenplot = pretty.time(thismonth,n=8)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for diurnal cycle...                                     #
   #---------------------------------------------------------------------------------------#
   thisday = seq(from=0,to=ndcycle,by=1) * 24 / ndcycle
   dcycplot = list()
   dcycplot$levels = c(0,4,8,12,16,20,24)
   dcycplot$n      = 7
   dcycplot$scale  = "hours"
   dcycplot$padj   = rep(0,times=dcycplot$n)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for soil profile layers...                               #
   #---------------------------------------------------------------------------------------#
   znice  = -pretty.log(-slz,n=8)
   znice  = sort(c(znice,slz[1],slz[nzg]))
   sel    = znice >= slz[1] & znice <= slz[nzg]
   znice  = znice[sel]
   zat    = -log(-znice)
   nznice = length(znice)
   znice  = sprintf("%.2f",znice)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define a suitable scale for monthly means...                                     #
   #---------------------------------------------------------------------------------------#
   montmont  = seq(from=1,to=12,by=1)
   montplot  = list()
   montplot$levels = montmont
   montplot$labels = capwords(mon2mmm(montmont))
   montplot$n      = 12
   montplot$scale  = "months"
   montplot$padj   = rep(0,times=dcycplot$n)
   #---------------------------------------------------------------------------------------#






   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #      Plotting section begins here...                                                  #
   #---------------------------------------------------------------------------------------#
   print ("    - Plotting figures...")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by PFT.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ntspft){
      thistspft   = tspft[[v]]
      vnam        = thistspft$vnam
      description = thistspft$desc
      unit        = thistspft$unit
      plotit      = thistspft$plt

      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && any(selpft)){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tspft",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"time series for all PFTs..."))

         #----- Load variable -------------------------------------------------------------#
         thisvar = get(vnam)

         #----- Loop over output formats. -------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if
            #---- Find the limit, and nudge it in case it is constant. --------------------#
            ylimit  = range(thisvar[,selpft],na.rm=TRUE)
            if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
               ylimit = c(-1,1)
            }else if(ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }#end if

            letitre = paste(description,lieu,sep=" - ")
            cols    = pftcols[selpft]
            legs    = pftnames[selpft]
            plot(x=thismonth,y=thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",xaxt="n",ylab=unit,cex.main=0.7)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){ 
               abline(v=whenplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
            for (n in 1:(npft+1)){
               if (selpft[n]){
                  lines(thismonth,thisvar[,n],type="l",col=pftcols[n],lwd=lwidth)
               }#end if
            }#end for
            legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Time series by LU.                                                               #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ntslu){
      thistslu    = tslu[[v]]
      vnam        = thistslu$vnam
      description = thistslu$desc
      unit        = thistslu$unit
      plotit      = thistslu$plt

      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit && any(sellu)){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tslu",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"time series for all LUs..."))

         #----- Load variable -------------------------------------------------------------#
         thisvar = get(vnam)

         #----- Loop over output formats. -------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if
            #---- Find the limit, and nudge it in case it is constant. --------------------#
            ylimit  = range(thisvar[,sellu],na.rm=TRUE)
            if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
               ylimit = c(-1,1)
            }else if(ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }#end if

            letitre = paste(description,lieu,sep=" - ")
            cols    = lucols[sellu]
            legs    = lunames[sellu]
            plot(thismonth,thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",ylab=unit,xaxt="n",cex.main=0.7)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){ 
               abline(v=whenplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
            for (n in 1:(nlu+1)){
               if (sellu[n]){
                  lines(thismonth,thisvar[,n],type="l",col=lucols[n],lwd=lwidth)
               }#end if
            }#end for
            legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot disturbance rate by disturbance transition.                                    #
   #---------------------------------------------------------------------------------------#
   if (tserdist && any(seldist)){
      print (paste("      + Disturbance rate time series for all disturbances..."))
      for (o in 1:nout){
         fichier = paste(outpref,"/disturb-",suffix,".",outform[o],sep="")
         if (outform[o] == "x11"){
            X11(width=size$width,height=size$height,pointsize=ptsz)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=size$width*depth,height=size$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=size$width,height=size$height
                      ,pointsize=ptsz,paper=paper)
         }#end if
         ylimit  = NULL
         for (jlu in 1:nlu){
            for (ilu in 1:nlu){
               if (seldist[ilu,jlu]) ylimit = range(c(ylimit,dist[,ilu,jlu]),na.rm=TRUE)
            }#end for
         }#end for
         if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
            ylimit = c(-1,1)
         }else if(ylimit[1] == ylimit[2] ){
            ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
            ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
         }#end if
         letitre = paste("Disturbance rates",lieu,sep=" - ")
         cols    = NULL
         legs    = NULL
         plot(thismonth,dist[,1,1],type="n",main=letitre,ylim=ylimit
             ,xlab="Time",ylab="[1/yr]",xaxt="n",cex.main=0.7)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){ 
               abline(v=whenplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
         n = 0
         for (jlu in 1:nlu){
            for (ilu in 1:nlu){
               n = n + 1
               if (seldist[ilu,jlu]){
                  cols = c(cols,distcols[n])
                  legs = c(legs,distnames[n])
                  lines(thismonth,dist[,ilu,jlu],type="l",col=distcols[n],lwd=lwidth)
               }#end if
            }#end for
         }#end for
         legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
      } #end for outform
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   print(paste("      * Plot some time series figures..."))
   for (hh in 1:ntser){

      #----- Retrieve variable information from the list. ---------------------------------#
      tsernow      = tser[[hh]]
      vnames       = tsernow$vnam  
      description  = tsernow$desc  
      lcolours     = tsernow$colour
      llwd         = tsernow$lwd
      ltype        = tsernow$type
      plog         = tsernow$plog
      prefix       = tsernow$prefix
      theme        = tsernow$theme 
      unit         = tsernow$unit  
      legpos       = tsernow$legpos
      plotit       = tsernow$plt   
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir = paste(outpref,"tseries",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"time series for several variables..."))


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         ylimit    = NULL
         for (l in 1:nlayers){
            thisvar = get(vnames[l])
            ylimit  = range(c(ylimit,thisvar),na.rm=TRUE)
         }#end for
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         print (paste("        > ",theme," time series ...",sep=""))

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",prefix,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if

            #----- Load variable ----------------------------------------------------------#
            thisvar = get(vnames[1])

            letitre = paste(theme," - ",lieu," \n"," Time series: ",theme,sep="")

            plot(x=thismonth,y=thisvar,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){ 
               abline(v=whenplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
            for (l in 1:nlayers){
               thisvar = get(vnames[l])
               points(x=thismonth,y=thisvar,col=lcolours[l]
                     ,lwd=llwd[l],type=ltype,pch=16,cex=0.8)
            }#end for
            legend(x=legpos,inset=0.05,legend=description,col=lcolours,lwd=llwd)
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the mean diurnal cycle.                                     #
   #---------------------------------------------------------------------------------------#
   print(paste("      * Plot some climatology of diurnal cycle..."))
   for (hh in 1:nclim){

      #----- Retrieve variable information from the list. ---------------------------------#
      climnow      = clim[[hh]]
      vnames       = climnow$vnam  
      description  = climnow$desc  
      lcolours     = climnow$colour
      llwd         = climnow$lwd
      ltype        = climnow$type
      plog         = climnow$plog
      prefix       = climnow$prefix
      theme        = climnow$theme 
      unit         = climnow$unit  
      legpos       = climnow$legpos
      plotit       = climnow$plt   
   
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"climdcyc",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = paste(outdir,prefix,sep="/")
         if (! file.exists(outtheme)) dir.create(outtheme)
         print (paste("      +",description,"diurnal cycle for several variables..."))


         #----- Define the number of layers. ----------------------------------------------#
         nlayers   = length(vnames)
         ylimit    = NULL
         for (l in 1:nlayers){
            thisvar = dcyc12mn[[vnames[l]]]
            ylimit  = range(c(ylimit,thisvar),na.rm=TRUE)
         }#end for
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in 1:12){
            cmon    = substring(100+pmon,2,3)
            namemon = mlist[pmon]

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            print (paste("        > ",theme," time series - ",namemon,"...",sep=""))

            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outtheme,"/",prefix,"-",cmon,"-",suffix,".",outform[o]
                              ,sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if

               #----- Load variable -------------------------------------------------------#
               thisvar = dcyc12mn[[vnames[1]]]
               thisvar = cbind(thisvar[,ndcycle],thisvar)

               letitre = paste(theme," - ",lieu,"\n"
                              ,"Mean diurnal cycle - ",namemon,sep="")

               plot(x=thisday,y=thisvar[pmon,],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=dcycplot$levels,labels=dcycplot$labels,padj=dcycplot$padj)
               if (plotgrid){ 
                  abline(v=dcycplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
               }#end if
               for (l in 1:nlayers){
                  thisvar = dcyc12mn[[vnames[l]]]
                  thisvar = cbind(thisvar[,ndcycle],thisvar)
                  points(x=thisday,y=thisvar[pmon,],col=lcolours[l]
                        ,lwd=llwd[l],type=ltype,pch=16,cex=0.8)
               }#end for
               legend(x=legpos,inset=0.05,legend=description,col=lcolours,lwd=llwd)
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end for pmon
      }#end if plotit
   }#end for ntser
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   print(paste("      * Comparisons of mean diurnal cycle (model vs. observations)..."))
   for (cc in 1:ncompdcyc){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compdcyc[[cc]]
      vname        = compnow$vnam  
      description  = compnow$desc  
      unit         = compnow$unit  
      plotsd       = compnow$plotsd
      lcolours     = compnow$colour
      errcolours   = compnow$errcol
      angle        = compnow$angle
      dens         = compnow$dens
      llwd         = compnow$lwd
      shwd         = compnow$shwd
      llwd         = compnow$lwd
      ltype        = compnow$type
      plog         = compnow$plog
      legpos       = compnow$legpos
      plotit       = compnow$plt

      #----- Check whether there are observations for this particular site. ---------------#
      if (iata == "mao"){
         obsnow = "obs.m34"
      }else if(iata == "stm"){
         obsnow = "obs.s67"
      }else if(iata == "rao"){
         obsnow = "obs.pdg"
      }else if(iata == "jpr"){
         obsnow = "obs.fns"
      }else if(iata == "btr"){
         obsnow = "obs.s77"
      }else{
         obsnow = paste("obs.",iata,sep="")
      }#end if
      plotit       = plotit && obsnow %in% ls()

      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisobs = get(obsnow)
         mnvar   = paste("qmean",vname,sep=".")
         sdvar   = paste("qsdev",vname,sep=".")
         obsmean = thisobs[[mnvar]]
         obssdev = thisobs[[sdvar]]
         #----- Append 1st hour after the last. -------------------------------------------#
         obsmean = cbind(obsmean,obsmean[,1])
         obssdev = cbind(obssdev,obssdev[,1])
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compdcyc",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         outtheme = paste(outdir,vname,sep="/")
         if (! file.exists(outtheme)) dir.create(outtheme)
         print (paste("      +",description,"comparison..."))
         #---------------------------------------------------------------------------------# 



         #----- Define the number of layers. ----------------------------------------------#
         thismean  = dcyc12mn[[vname]]
         thissdev  = dcyc12sd[[vname]]
         #---------------------------------------------------------------------------------# 



         #---------------------------------------------------------------------------------# 
         #    Some variables have no standard deviation in the model.  Make them 0 if this #
         # is the case.                                                                    #
         #---------------------------------------------------------------------------------# 
         if (length(thissdev) == 0){
            thissdev = 0. * thismean
         }#end if
         #---------------------------------------------------------------------------------# 


         #----- Append the last hour before the first one. --------------------------------#
         thismean  = cbind(thismean[,ndcycle],thismean)
         thissdev  = cbind(thissdev[,ndcycle],thissdev)
         #---------------------------------------------------------------------------------# 


         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = range(c(thismean + thissdev ,thismean - thissdev
                               ,obsmean  + obssdev  ,obsmean  - obssdev    ),na.rm=TRUE)
         }else{
            ylimit    = range(c(thismean,obsmean),na.rm=TRUE)
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all months.                                                      #
         #---------------------------------------------------------------------------------#
         for (pmon in 1:12){
            cmon    = substring(100+pmon,2,3)
            namemon = mlist[pmon]

            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            print (paste("        > ",description," time series - ",namemon,"...",sep=""))

            #----- Loop over formats. -----------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outtheme,"/",vname,"-",cmon,".",outform[o]
                              ,sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if

               #----- Load variable -------------------------------------------------------#
               letitre = paste(description," - ",lieu,"\n"
                              ,"Mean diurnal cycle - ",namemon,sep="")
               plot(x=thisday,y=thismean[pmon,],type="n",main=letitre,xlab="Time"
                   ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                   ,cex.main=cex.main)
               axis(side=1,at=dcycplot$levels,labels=dcycplot$labels,padj=dcycplot$padj)
               if (plotgrid){ 
                  abline(v=dcycplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
               }#end if
               if (plotsd){
                  err.x = c(thisday,rev(thisday),NA,thisday,rev(thisday))
                  err.y = c(thismean[pmon,] + thissdev[pmon,]
                           ,rev(thismean[pmon,]) - rev(thissdev[pmon,])
                           ,NA
                           ,obsmean[pmon,]      + obssdev[pmon,]
                           ,rev(obsmean[pmon,]) - rev(obssdev[pmon,]))
                  polygon(x=err.x,y=err.y,col=errcolours,angle=angle,density=dens
                         ,lty="solid",lwd=shwd)
               }#end if
               points(x=thisday,y=thismean[pmon,],col=lcolours[1]
                     ,lwd=llwd[1],type=ltype,pch=16,cex=1.0)
               points(x=thisday,y=obsmean[pmon,],col=lcolours[2]
                     ,lwd=llwd[2],type=ltype,pch=16,cex=1.0)
               if (plotsd){
                  legend(x=legpos,inset=0.05,legend=c("Model","Observation")
                        ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                        ,bg="white",title="Shaded areas = 1 SD",cex=1.0,pch=16)
               }else{
                  legend(x=legpos,inset=0.05,legend=c("Model","Observation")
                        ,col=lcolours,lwd=llwd,cex=1.0,pch=16)
               }#end if
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end for pmon
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the comparison between observations and model.                                 #
   #---------------------------------------------------------------------------------------#
   print(paste("    + Comparisons of monthly means (model vs. observations)..."))
   for (cc in 1:ncompmmean){

      #----- Retrieve variable information from the list. ---------------------------------#
      compnow      = compmmean[[cc]]
      vname        = compnow$vnam  
      description  = compnow$desc  
      unit         = compnow$unit  
      plotsd       = compnow$plotsd
      lcolours     = compnow$colour
      errcolours   = compnow$errcol
      angle        = compnow$angle
      dens         = compnow$dens
      llwd         = compnow$lwd
      shwd         = compnow$shwd
      llwd         = compnow$lwd
      ltype        = compnow$type
      plog         = compnow$plog
      legpos       = compnow$legpos
      plotit       = compnow$plt

      #----- Check whether there are observations for this particular site. ---------------#
      if (iata == "mao"){
         obsnow = "obs.m34"
      }else if(iata == "stm"){
         obsnow = "obs.s67"
      }else if(iata == "rao"){
         obsnow = "obs.pdg"
      }else if(iata == "jpr"){
         obsnow = "obs.fns"
      }else if(iata == "btr"){
         obsnow = "obs.s77"
      }else{
         obsnow = paste("obs.",iata,sep="")
      }#end if

      plotit       = plotit && obsnow %in% ls()





      if (plotit){
         #---------------------------------------------------------------------------------#
         #    Copy the observations to a scratch variable.                                 #
         #---------------------------------------------------------------------------------#
         thisobs = get(obsnow)
         mnvar   = paste("mmean",vname,sep=".")
         sdvar   = paste("msdev",vname,sep=".")
         obsmean = thisobs[[mnvar]]
         obssdev = thisobs[[sdvar]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Check whether the time series directory exists.  If not, create it.          #
         #---------------------------------------------------------------------------------#
         outdir   = paste(outpref,"compmmean",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      - ",description,"comparison..."))
         #---------------------------------------------------------------------------------#



         #----- Define the number of layers. ----------------------------------------------#
         thismean  = mont12mn[[vname]]
         thissdev  = mont12sd[[vname]]
         #---------------------------------------------------------------------------------# 



         #---------------------------------------------------------------------------------# 
         #    Some variables have no standard deviation in the model.  Make them 0 if this #
         # is the case.                                                                    #
         #---------------------------------------------------------------------------------# 
         if (length(thissdev) == 0){
            thissdev = 0. * thismean
         }#end if
         #---------------------------------------------------------------------------------# 



         #----- Find the plot range. ------------------------------------------------------#
         if (plotsd){
            ylimit    = range(c(thismean + thissdev ,thismean - thissdev
                               ,obsmean  + obssdev  ,obsmean  - obssdev    ),na.rm=TRUE)
         }else{
            ylimit    = range(c(thismean,obsmean),na.rm=TRUE)
         }#end if
         #----- Expand the upper range in so the legend doesn't hide things. --------------#
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vname,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if

            #----- Load variable ----------------------------------------------------------#
            letitre = paste(description," - ",lieu,"\n","Monthly mean",sep="")
            plot(x=montmont,y=thismean,type="n",main=letitre,xlab="Time"
                ,ylim=ylimit,ylab=paste("[",unit,"]",sep=""),log=plog,xaxt="n"
                ,cex.main=cex.main)
            axis(side=1,at=montplot$levels,labels=montplot$labels,padj=montplot$padj)
            if (plotgrid){ 
               abline(v=montplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
            if (plotsd){
               err.x = c(montmont,rev(montmont),NA,montmont,rev(montmont))
               err.y = c(thismean + thissdev,rev(thismean) - rev(thissdev),NA
                        ,obsmean  + obssdev ,rev(obsmean ) - rev(obssdev )   )
               polygon(x=err.x,y=err.y,col=errcolours,angle=angle,density=dens
                      ,lty="solid",lwd=shwd)
            }#end if
            points(x=montmont,y=thismean,col=lcolours[1],lwd=llwd[1],type=ltype
                  ,pch=16,cex=1.0)
            points(x=montmont,y=obsmean ,col=lcolours[2],lwd=llwd[2],type=ltype
                  ,pch=16,cex=1.0)
            if (plotsd){
               legend(x=legpos,inset=0.05,legend=c("Model","Observation")
                     ,fill=errcolours,angle=angle,density=dens,lwd=llwd,col=lcolours
                     ,bg="white",title="Shaded areas = 1 SD",cex=1.0,pch=16)
            }else{
               legend(x=legpos,inset=0.05,legend=c("Model","Observation")
                     ,col=lcolours,lwd=llwd,cex=1.0,pch=16)
            }#end if
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for ncompare
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the climatology of the soil properties.                                        #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nsoilclim){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisclim    = soilclim[[v]]
      vnam        = thisclim$vnam
      description = thisclim$desc
      unit        = thisclim$unit
      vcscheme    = thisclim$csch
      pnlog       = thisclim$pnlog
      plotit      = thisclim$plt

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"soilclim",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      + Climatology profile of ",description,"..."))

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis  = sort(unique(monnum))
         soilaxis = slz
         nmon     = length(monaxis)
         nsoil    = nzg

         #----- Save the meaningful months and years. -------------------------------------#
         monat   = 1:12
         monlab  = c("J","F","M","A","M","J","J","A","S","O","N","D")

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = mont12mn[[vnam]]

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         january  = vararr[1,]
         january  = c(january,january[nzg],january[nzg])

         december = vararr[12,]
         december = c(december[1],december[1],december)

         #----- Bind first and last year to the array, to make the edges buffered. ---------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nzg])
         varbuff  = rbind(december,varbuff,january)

         #----------------------------------------------------------------------------------#
         #   Expand the month and year axes.  Make the -------------------------------------------#
         monaxis  = c(min(monaxis)-1,monaxis,max(monaxis)+1)
         soilaxis = -log(-1.0 * c( slz[1]*(slz[1]/slz[2])
                                 , soilaxis
                                 , slz[nzg]*(slz[nzg]/slz[nzg-1]) ))

         if (pnlog){
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty.log(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)
         }else{
            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)
         }#end if

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if

            letitre = paste(description," - ",lieu,sep="")
            sombreado(x=monaxis,y=soilaxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Soil depth [m]"
                                      ,cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,key.log=pnlog
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=zat,labels=znice)
                                 if (hovgrid){
                                    abline(h=zat,v=monat,col="lightgray",lty="dotted")
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the Hovmoller looking diagrams showing months and years.                       #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nhov){

      #----- Retrieve variable information from the list. ---------------------------------#
      thishovdi   = hovdi[[v]]
      vnam        = thishovdi$vnam
      description = thishovdi$desc
      unit        = thishovdi$unit
      vcscheme    = thishovdi$csch
      plotit      = thishovdi$plt

      #------------------------------------------------------------------------------------#
      #     Find the first and the last full years.  These will be the actual first and    #
      # last year only if the years are complete, otherwise the first and the last year    #
      # will be taken out.                                                                 #
      #------------------------------------------------------------------------------------#
      if (monthbeg == 1){
         yearaa = yeara
      }else{
         yearaa = yeara + 1
      }# end if
      if (meszz == 12){
         yearzz = yearz
      }else{
         yearzz = yearz - 1
      }#end if
      sel      = myear >= yearaa & myear <= yearzz
      twoyears = sum(sel) >= 24

      if (plotit && twoyears){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"hovmoller",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"Hovmoller time series ..."))

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = get(vnam)

         #----- Find the number of rows and columns, and the axes. ------------------------#
         monaxis = sort(unique(monnum[sel]))
         yraxis  = sort(unique(myear[sel]))
         nmon    = length(monaxis)
         nyear   = length(yraxis)

         #----- Save the meaningful months and years. -------------------------------------#
         monat   = 1:12
         monlab  = c("J","F","M","A","M","J","J","A","S","O","N","D")
         yrat    = pretty(yraxis)

         #----- Convert the vector data into an array. ------------------------------------#
         vararr  = array(thisvar[sel],c(nmon,nyear))

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         january  = vararr[1,]
         january  = c(january,january[nyear],january[nyear])

         december = vararr[12,]
         december = c(december[1],december[1],december)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = cbind(vararr[,1],vararr,vararr[,nyear])
         varbuff  = rbind(december,varbuff,january)

         #----- Expand the month and year axes. -------------------------------------------#
         monaxis = c(min(monaxis)-1,monaxis,max(monaxis)+1)
         yraxis  = c(min(yraxis)-1,yraxis,max(yraxis)+1)

         vrange  = range(varbuff,na.rm=TRUE)
         vlevels = pretty(x=vrange,n=ncolshov)
         vnlev   = length(vlevels)

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if

            letitre = paste(description," - ",lieu,sep="")
            sombreado(x=monaxis,y=yraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,xlab="Month",ylab="Year",cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                 axis(side=2,at=yrat)
                                 if (hovgrid){
                                    for (yl in yrat){
                                       abline(h=yl,col="lightgray",lty="dotted")
                                    } #end for yl
                                    for (ml in monat){
                                       abline(v=ml,col="lightgray",lty="dotted")
                                    } #end for ml
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #   Plot the Hovmoller looking diagrams showing time of day and time.                   #
   #---------------------------------------------------------------------------------------#

   for (v in 1:nhdcyc){

      #----- Retrieve variable information from the list. ---------------------------------#
      thishdcyc   = hdcyc[[v]]
      vnam        = thishdcyc$vnam
      description = thishdcyc$desc
      unit        = thishdcyc$unit
      vcscheme    = thishdcyc$csch
      plotit      = thishdcyc$plt

      if (plotit){

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"hovdcycle",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"time series of diurnal cycle..."))

         #----- Load this variable into "thisvar". ----------------------------------------#
         vararr   = dcycmean[[vnam]]

         #----- Copy Decembers ans Januaries to make the edges buffered. ------------------#
         firsthr  = vararr[,1]
         firsthr  = c(firsthr,firsthr[totmon],firsthr[totmon])

         lasthr   = vararr[,ndcycle]
         lasthr   = c(lasthr[1],lasthr[1],lasthr)

         #----- Bind first and last year to the array, to make the edges buffered. --------#
         varbuff  = rbind(vararr[1,],vararr,vararr[totmon,])
         varbuff  = cbind(lasthr,varbuff,firsthr)

         #----- Expand the month and year axes. -------------------------------------------#
         hraxis    = seq(from=0,to=ndcycle+1,by=1) * 24 / ndcycle
         dwhen     = thismonth[2]-thismonth[1]
         whenaxis  = c(thismonth[1]-dwhen,thismonth,thismonth[totmon]+dwhen)
         hdcycplot = pretty.time(whenaxis,n=8)

         vrange  = range(varbuff,na.rm=TRUE)
         vlevels = pretty(x=vrange,n=ncolshov)
         vnlev   = length(vlevels)

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if

            letitre = paste("Mean diurnal cycle \n ",description," - ",lieu,sep="")
            sombreado(x=whenaxis,y=hraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                     ,color.palette=get(vcscheme)
                     ,plot.title=title(main=letitre,ylab="Time of day (GMT)"
                                      ,xlab="Time",cex.main=0.7)
                     ,key.title=title(main=unit,cex.main=0.8)
                     ,plot.axes={axis(side=1,at=hdcycplot$level,labels=hdcycplot$labels)
                                 axis(side=2,at=dcycplot$levels,labels=dcycplot$labels)
                                 if (hovgrid){
                                    abline(v=hdcycplot$levels,h=dcycplot$levels
                                          ,col="lightgray",lty="dotted")
                                 }#end if hovgrid
                                }#end plot.axes
                     )

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Plot the monthly boxplots.                                                          #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nbox){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisbplot   = bplot[[v]]
      vnam        = thisbplot$vnam
      description = thisbplot$desc
      unit        = thisbplot$unit
      plotit      = thisbplot$plt

      if (plotit){
         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"boxplot",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",description,"box plot..."))

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = get(vnam)

         for (o in 1:nout){
            fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
            if (outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if
            ylimit  = range(thisvar, na.rm=TRUE)
            if (ylimit[1] == ylimit[2] && ylimit[1] == 0){
               ylimit = c(-1,1)
            }else if(ylimit[1] == ylimit[2] ){
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * ylnudge)
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * ylnudge)
            }#end if
            letitre = paste(description,lieu,sep=" - ")
            plot(mmonth,thisvar,main=letitre,ylim=ylimit,cex.main=0.7
                ,xlab="Time",ylab=paste("[",unit,"]",sep=""))

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if
   }#end for nbox
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #    Plot the 3-D size and age structure of light level.                                #
   #---------------------------------------------------------------------------------------#
   for (v in 1:npsas){
      #----- Retrieve variable information from the list. ---------------------------------#
      thissas     = sasplot[[v]]
      vnam        = thissas$vnam
      description = thissas$desc
      unit        = thissas$unit
      plotit      = thissas$plt

      #----- If this variable is to be plotted, then go through this if block. ------------#
      if (plotit){

         print (paste("      +",description,"size and age structure plot..."))

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         sasdir = paste(outpref,"sas",sep="/")
         if (! file.exists(sasdir)) dir.create(sasdir)
         outdir = paste(sasdir,vnam,sep="/")
         if (! file.exists(outdir)) dir.create(outdir)

         #----- Load this list into "thislist". -------------------------------------------#
         varco =  get(vnam)


         for (ww in names(ageco)){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth   = substring(ww,7,8)
            thisyear = substring(ww,2,5)
            mm       = as.numeric(cmonth)
            yy       = as.numeric(thisyear)

            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            ageww   = ageco[[ww]]
            dbhww   = dbhco[[ww]]
            pftww   = pftco[[ww]]
            varww   = varco[[ww]]
            popww   = nplantco[[ww]] * areaco[[ww]]

            #------------------------------------------------------------------------------#
            #     We only plot the SAS figures when the polygon is not an absolute desert. #
            #------------------------------------------------------------------------------#
            if (any (! is.na(varww))){
               #---------------------------------------------------------------------------#
               #      Find the range.  If the user wants the range to be fixed, then use   #
               # the global range, otherwise, simply use the range for this year.          #
               #---------------------------------------------------------------------------#
               if (sasfixlimits){
                  xlimit  = range(ageco            , na.rm=TRUE)
                  ylimit  = range(dbhco            , na.rm=TRUE)
                  zlimit  = range(varco            , na.rm=TRUE)
                  popmin  = min  (nplantco * areaco, na.rm=TRUE)
                  popmax  = max  (nplantco * areaco, na.rm=TRUE)
               }else{
                  xlimit  = range(ageww  ,na.rm=TRUE)
                  ylimit  = range(dbhww  ,na.rm=TRUE)
                  zlimit  = range(varww  ,na.rm=TRUE)
                  popmin  = min  (popww  ,na.rm=TRUE)
                  popmax  = max  (popww  ,na.rm=TRUE)
               }#end if

               #----- Define the scale-dependent population size. -------------------------#
               cexww = cexmin + (cexmax - cexmin) * log(popww/popmin) / log(popmax/popmin)

               #----- Define the floor location. ------------------------------------------#
               if (zlimit[1] == zlimit[2]){
                  if (zlimit[1] == 0){
                     zlimit = c(-1.,1.)
                  }else{
                     zlimit = sort(c(0.9,1.1)*zlimit[1])
                  }#end if
               }#end if
               if ((zlimit[1] > 0) != (zlimit[2] > 0)){
                  floor3d = 0.
               }else if (zlimit[1] > 0){
                  floor3d = zlimit[1]
               }else{
                  floor3d = zlimit[2]
               }#end if

               #----- Define the grid information for the 3-D plot. -----------------------#
               ageaxis   = pretty(xlimit,n=20)
               dbhaxis   = pretty(ylimit,n=20)
               xlimit    = range(ageaxis)
               ylimit    = range(dbhaxis)
               flooraxis = matrix(floor3d,nrow=length(ageaxis),ncol=length(dbhaxis))

               #----- Expand the lines to make the lollipops. -----------------------------#
               ncohnow  = length(varww)
               ageww    = rep(ageww,each=3)
               dbhww    = rep(dbhww,each=3)
               pftww    = rep(pftww,each=3)
               varww    = as.vector(rbind(rep(floor3d,times=ncohnow)
                                         ,varco[[ww]]
                                         ,rep(NA,times=ncohnow)))
               pchww    = rep(c(NA,16,NA),times=ncohnow)
               cexww    = rep(cexww,each=3)
               colww    = pftcols[pftww]

               pftin   = sort(unique(pftco[[ww]]))
               colleg  = pftcols[pftin]
               pftleg  = pftnames[pftin]


               #----- Loop over output formats. -------------------------------------------#
               for (o in 1:nout){
                  fichier = paste(outdir,"/",vnam,"-",thisyear,"-",cmonth,"-",suffix
                                            ,".",outform[o],sep="")
                  if (outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=paper)
                  }#end if

                  stcol   = pftcols[pftww]
                  letitre = paste(description," - ",lieu,
                                  "\n Time :",mlist[mm],"/",thisyear,sep=" ")
                  lezlab  = paste(description," [",unit,"]",sep="")

                  #----- First plot: the box. ---------------------------------------------#
                  pout = persp(x=ageaxis,y=dbhaxis,z=flooraxis,xlim=xlimit,ylim=ylimit
                              ,zlim=zlimit,theta=theta,phi=phi,col=gcol,expand=expz
                              ,ticktype="detailed",border=NA,xlab="Gap age [yr]"
                              ,ylab="DBH [cm]",zlab=lezlab,shade=shade,ltheta=ltheta
                              ,main=letitre,cex.main=0.7)
                  #----- Second plot, the actual data (aka my lollipop trees). ------------#
                  lines (trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="l"
                        ,col="gray29",lwd=2)
                  points(trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="p"
                        ,pch=pchww,col=colww,cex=cexww)
                  legend(x="bottomright",inset=0.01,legend=pftleg,fill=colleg
                        ,ncol=1,bg="white",cex=0.9)


                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
               } #end for outform
            }#end if is.na(varww)
         }#end for nameco
      } #end if
   }#end for npsas
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the filled contour plots as function of time and PFT.                          #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nfcpft){

      #----- Retrieve variable information from the list. ---------------------------------#
      thisfc      = fcpft[[v]]
      vnam        = thisfc$vnam
      description = thisfc$desc
      unit        = thisfc$unit
      vcscheme    = thisfc$csch
      vclass      = thisfc$cls
      plotit      = thisfc$plt

      if (plotit && any(selpft)){
         #---------------------------------------------------------------------------------#
         #      Define which class we are going to plot.                                   #
         #---------------------------------------------------------------------------------#
         if (vclass == "age"){
            thisclass = classage
            expclass  = c(classage[1]-dage,classage,classage[nage]+dage)
            nclass    = nage
            atclass   = pretty(classage)
            leylab    = "Age class"
         }else{ # if (vclass == "dbh"){
            thisclass = classdbh
            expclass  = c(classdbh[1]-ddbh,classdbh,classdbh[ndbh]+ddbh)
            nclass    = ndbh
            atclass   = pretty(classdbh)
            leylab    = "DBH class"
         }#end if vclass

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir  =  paste(outpref,"fcpft",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      + ",description," time series by ",leylab,"...",sep=""))

         #----- Load this variable into "thisvar". ----------------------------------------#
         thisvar = get(vnam)

         #---------------------------------------------------------------------------------#
         #     Define which PFTs to run.  Also, find the overall range amongst all PFTs,   #
         # so all plants are in the same scale.                                            #
         #---------------------------------------------------------------------------------#
         selpftonly = selpft[1:npft]
         pftrun  = seq(from=1,to=npft,by=1)[selpftonly]
         vrange  = range(thisvar[,,selpftonly],na.rm=TRUE)
         vlevels = pretty(x=vrange,n=ncolsfc)
         vnlev   = length(vlevels)


         for (p in pftrun){
            if (p < 10){
               cpcp = paste("0",p,sep="")
            }else{
               cpcp = as.character(p)
            } #end if
            for (o in 1:nout){
               fichier = paste(outdir,"/",vnam,"-pft",cpcp,"-",suffix,".",outform[o]
                              ,sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if

               #----- Copy this PFT to a scratch array. -----------------------------------#
               var00 = as.matrix(thisvar[,,p])
               
               #----- Expand the axis so all meaningful information is in the middle. -----#
               before   = thismonth[1]-(thismonth[2]-thismonth[1])
               after    = thismonth[tlast]+(thismonth[tlast]-thismonth[tlast-1])
               expmonth = c(before,thismonth,after)
               expmonth = chron(expmonth,out.format=c(dates="day-mon-yr",times=NULL))
               zero     = rep(0,times=tlast)
               expvar   = cbind(zero,var00,zero)
               expvar   = rbind(expvar[1,],expvar,expvar[tlast,])

               letitre = paste(description,".  PFT:",p," - ",lieu,sep="")
               sombreado(x=expmonth,y=expclass,z=expvar,levels=vlevels,nlevels=vnlev
                        ,color.palette=get(vcscheme)
                        ,plot.title=title(main=letitre,xlab="Time",ylab=leylab,cex.main=0.7)
                        ,key.title=title(main=unit,cex.main=0.8)
                        ,plot.axes={axis(side=1,at=whenplot$levels,labels=whenplot$labels
                                        ,padj=whenplot$padj)
                                    axis(side=2,at=atclass)
                                    if (fcgrid){ 
                                       abline(v=whenplot$levels,h=atclass,col="lightgray"
                                             ,lty="solid")
                                    }#end if
                                   }#end plot.axes
                        )

               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end for pftrun
      }#end if plotit
   }#end for nfcpft
   #---------------------------------------------------------------------------------------#
}#end for places

#q("no")
