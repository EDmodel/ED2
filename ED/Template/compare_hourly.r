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
reload.hour = TRUE
rdata.path  = file.path(here,"RData_often")
rdata.hour  = file.path(rdata.path,"hourly_ed22.RData")
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Site settings:                                                                       #
# eort      -- first letter ("e" or "t")                                                   #
# sites     -- site codes ("IATA")                                                         #
# sites$pch -- site symbols                                                                #
#------------------------------------------------------------------------------------------#
eort           = "t"
sites          = list( list( iata = "gyf", desc = "Paracou"         , pch =  2)
                     , list( iata = "s67", desc = "Santarem (km 67)", pch =  5)
                     , list( iata = "s83", desc = "Santarem (km 83)", pch =  9)
                     , list( iata = "pdg", desc = "Pe-de-Gigante"   , pch = 13)
                     , list( iata = "rja", desc = "Rebio Jaru"      , pch =  1)
                     , list( iata = "m34", desc = "Manaus (K34)"    , pch =  6)
#                    , list( iata = "pnz", desc = "Petrolina"       , pch =  4)
#                    , list( iata = "ban", desc = "Bananal"         , pch =  8)
#                    , list( iata = "cax", desc = "Caxiuana"        , pch =  0)
                     )#end list
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.struct     = list( name     = c("ble_iage30_pft02","ble_iage30_pft05"
                                   ,"sas_iage01_pft02","sas_iage01_pft05"
                                   ,"sas_iage30_pft02","sas_iage30_pft05"
                                   )#end c
                     , desc     = c("Size 01 + Age 01 + PFT 02"
                                   ,"Size 01 + Age 01 + PFT 05"
                                   ,"Size 80 + Age 01 + PFT 02"
                                   ,"Size 80 + Age 01 + PFT 05"
                                   ,"Size 80 + Age 30 + PFT 02"
                                   ,"Size 80 + Age 30 + PFT 05"
                                   )#end c
                     , verbose  = c("2 PFTs"             ,"5 PFTs"
                                   ,"Size + 2 PFTs"      ,"Size + 5 PFTs"
                                   ,"Size + Age + 2 PFTs","Size + Age + 5 PFTs"
                                   )#end c
                     , colour   = c("#6600FF","#CC00FF"
                                   ,"#FF3300","#FF6633"
                                   ,"#009900","#00FF66"
                                   )#end c
                     , fgcol    = c("#330080","#660080"
                                   ,"#801A00","#80331A"
                                   ,"#004D00","#008033"
                                   )#end c
                     )#end list
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
wpaper         = "legal"               # Wide paper size, to define the plot shape
ptsz           = 16                    # Font size.
lwidth         = 2.5                   # Line width
plotgrid       = TRUE                  # Should I plot the grid in the background? 

legwhere       = "topleft"             # Where should I place the legend?
inset          = 0.01                  # Inset between legend and edge of plot region.
fracexp        = 0.40                  # Expansion factor for y axis (to fit legend)
cex.main       = 0.8                   # Scale coefficient for the title

st.cex.min     = 0.7                   # Minimum and maximum sizes for points in the 
st.cex.max     = 2.0                   #     Skill and Taylor diagrams
st.lwd.min     = 1.3                   # Minimum and maximum sizes for points in the 
st.lwd.max     = 3.0                   #     Skill and Taylor diagrams


light.method   = "nls"                 # Which method to apply (nls or optim)
light.skew     = FALSE                 # Use skew-normal distribution on residuals?
light.n.boot   = 1000                  # # of bootstrap iterations (light response fit)
ftnight.n.boot = 50                    # # of bootstrap iterations (fn mean)

n.quant        = 1024                  # # of quantiles to produce the density function.
                                       #    We strongly advise to choose a number that is
                                       #    a power of two, especially when using EDF 
                                       #    (otherwise distributions will be interpolated).
nhour.min      = 16                    # Minimum number of hours to use the data.

#------------------------------------------------------------------------------------------#
#      Switch controls to plot only the needed ones.                                       #
#------------------------------------------------------------------------------------------#
plot.light         = c(FALSE,TRUE)[2]
plot.ts.ftnight    = c(FALSE,TRUE)[1]
plot.bp.diel       = c(FALSE,TRUE)[1]
plot.qq.dmean      = c(FALSE,TRUE)[1]
plot.density.dmean = c(FALSE,TRUE)[1]
plot.spider        = c(FALSE,TRUE)[1]
plot.skill.taylor  = c(FALSE,TRUE)[1]
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
compvar[[ n]] = list( vnam       = "ustar"
                    , symbol     = expression(u^symbol("\052"))
                    , desc       = "Friction velocity"
                    , unit       = untab$mos
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(purple.bg,purple.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "cflxca"
                    , symbol     = expression(F(CO[2]))
                    , desc       = "Carbon dioxide flux"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(green.bg,green.fg)
                    , leg.corner = "bottomright"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "cflxst"
                    , symbol     = expression(S(CO[2]))
                    , desc       = "Carbon dioxide storage"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(orange.bg,orange.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "nee"
                    , symbol     = expression(NEE)
                    , desc       = "Net ecosystem exchange"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(green.bg,green.fg)
                    , leg.corner = "bottomright"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "nep"
                    , symbol     = expression(NEP)
                    , desc       = "Net ecosystem productivity"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(olive.bg,olive.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "reco"
                    , symbol     = expression(R[Eco])
                    , desc       = "Ecosystem respiration"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(yellow.bg,yellow.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "gpp"
                    , symbol     = expression(GPP)
                    , desc       = "Gross primary productivity"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(green.bg,green.fg)
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "parup"
                    , symbol     = expression(PAR^symbol("\335"))
                    , desc       = "Outgoing PAR"
                    , unit       = untab$umolom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(olive.bg,olive.fg)
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "rshortup"
                    , symbol     = expression(SW^symbol("\335"))
                    , desc       = "Outgoing shortwave radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(indigo.bg,indigo.fg)
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
#n = n + 1
#compvar[[ n]] = list( vnam       = "rnet"
#                    , symbol     = expression(R[Net])
#                    , desc       = "Net radiation"
#                    , unit       = untab$wom2
#                    , col.obser  = c(grey.bg,grey.fg)
#                    , col.model  = c(sky.bg,sky.fg)
#                    , leg.corner = "topleft"
#                    , sunvar     = FALSE
#                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "rlongup"
                    , symbol     = expression(LW^symbol("\335"))
                    , desc       = "Outgoing longwave radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(red.bg,red.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "hflxca"
                    , symbol     = expression(F(theta))
                    , desc       = "Sensible heat flux"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(orange.bg,orange.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
n = n + 1
compvar[[ n]] = list( vnam       = "wflxca"
                    , symbol     = expression(F(H[2]*O))
                    , desc       = "Water vapour flux"
                    , unit       = untab$kgwom2oday
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(blue.bg,blue.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Input variables.                                                                     #
#------------------------------------------------------------------------------------------#
control       = list()
control[[ 1]] = list( vnam       = "rshort"
                    , desc       = "Incoming shortwave radiation"
                    , unit       = untab$wom2
                    )#end list
control[[ 2]] = list( vnam       = "rlong"
                    , desc       = "Incoming longwave radiation"
                    , unit       = untab$wom2
                    )#end list
control[[ 3]] = list( vnam       = "atm.prss"
                    , desc       = "Air pressure"
                    , unit       = untab$hpa
                    )#end list
control[[ 4]] = list( vnam       = "atm.temp"
                    , desc       = "Air temperature"
                    , unit       = untab$degC
                    )#end list
control[[ 5]] = list( vnam       = "atm.shv"
                    , desc       = "Air specific humidity"
                    , unit       = untab$gwokg
                    )#end list
control[[ 6]] = list( vnam       = "atm.vels"
                    , desc       = "Wind speed"
                    , unit       = untab$mos
                    )#end list
control[[ 7]] = list( vnam       = "rain"
                    , desc       = "Precipitation rate"
                    , unit       = untab$kgwom2oday
                    )#end list
control[[ 8]] = list( vnam       = "bsa"
                    , desc       = "Basal area"
                    , unit       = untab$cm2om2
                    )#end list
control[[ 9]] = list( vnam       = "wdens"
                    , desc       = "Mean wood density"
                    , unit       = untab$kgom3
                    )#end list
control[[10]] = list( vnam       = "global"
                    , desc       = "Global index"
                    , unit       = untab$empty
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
#      List the keys for all dimensions.                                                   #
#------------------------------------------------------------------------------------------#
vnum        = which(sapply(X=sites[[1]],FUN=is.numeric))
vlog        = which(sapply(X=sites[[1]],FUN=is.logical))
sites       = data.frame( apply(X=sapply(X=sites,FUN=c),MARGIN=1,FUN=unlist)
                        , stringsAsFactors = FALSE
                        )#end data.frame
for (vn in vnum) sites[[vn]] = as.numeric(sites[[vn]])
for (vl in vlog) sites[[vl]] = as.logical(sites[[vl]])

control.key = apply(X = sapply(X=control,FUN=c),MARGIN=1,FUN=unlist)[,"vnam"]
compvar.key = apply(X = sapply(X=compvar,FUN=c),MARGIN=1,FUN=unlist)$vnam
compvar.sym = apply(X = sapply(X=compvar,FUN=c),MARGIN=1,FUN=unlist)$symbol
good.key    = apply(X = sapply(X=good   ,FUN=c),MARGIN=1,FUN=unlist)[,"vnam"]
season.key  = season.list
diel.key    = c("night","rise.set","day","all.hrs","dmean")
diel.desc   = c("Nighttime","Sun Rise/Set","Daytime","All hours","Daily mean")
moment.key  = c("mean","variance","skewness","kurtosis")
moment.desc = c("Mean","Variance","Skewness","Kurtosis")
hour.num    = seq(from=0,to=23,by=3)
hour.key    = c( paste( sprintf("%2.2i",(hour.num - 1) %% 24)
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
nsites   = length(sites$iata )
nsimul   = length(simul.key  )
ncompvar = length(compvar.key)
ncontrol = length(control.key)
ngood    = length(good.key   )
nseason  = length(season.key )
ndiel    = length(diel.key   )
nhour    = length(hour.key   )
nmoment  = length(moment.key )
#------------------------------------------------------------------------------------------#

season.suffix = paste(sprintf("%2.2i",sequence(nseason)),tolower(season.key),sep="-")


#----- Load observations. -----------------------------------------------------------------#
obser.file = paste(srcdir,"LBA_MIP.nogapfill.RData",sep="/")
load(file=obser.file)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size  = plotsize(proje=FALSE,paper=paper )
wsize = plotsize(proje=FALSE,paper=wpaper)
#------------------------------------------------------------------------------------------#



#----- Find the best set up for plotting all seasons in the same plot. --------------------#
lo.box   = pretty.box(n=nseason-1)
lo.simul = pretty.box(n=nsimul,byrow=FALSE)
lo.site  = pretty.box(n=nsites)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
out = list()
for (o in 1:nout){
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
   o.spider                = list()
   o.spider$main           = file.path(o.form$main,"spider")
   if (is.figure && ! file.exists(o.spider$main)) dir.create(o.spider$main)
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
   o.ts.ftnight           = list()
   o.ts.ftnight$main      = file.path(o.form$main,"ts_ftnight")
   if (is.figure && ! file.exists(o.ts.ftnight$main)) dir.create(o.ts.ftnight$main)
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam

      o.compvar        = file.path(o.ts.ftnight$main,this.vnam)
      if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
      #------------------------------------------------------------------------------------#
      o.ts.ftnight[[this.vnam]] = o.compvar
   }#end for (v in sequence(ncompvar))
   o.form$ts.ftnight = o.ts.ftnight
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the box plots of the mean diel.                                  #
   #---------------------------------------------------------------------------------------#
   o.bp.diel           = list()
   o.bp.diel$main      = file.path(o.form$main,"bp_diel")
   if (is.figure && ! file.exists(o.bp.diel$main)) dir.create(o.bp.diel$main)
   for (v in sequence(ncompvar)){
      this.compvar     = compvar[[v]]
      this.vnam        = this.compvar$vnam

      o.compvar        = file.path(o.bp.diel$main,this.vnam)
      if (is.figure && ! file.exists(o.compvar)) dir.create(o.compvar)
      #------------------------------------------------------------------------------------#
      o.bp.diel[[this.vnam]] = o.compvar
   }#end for (v in sequence(ncompvar))
   o.form$bp.diel = o.bp.diel
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the Q-Q plots of daily means.                                    #
   #---------------------------------------------------------------------------------------#
   o.qq.dmean       = list()
   o.qq.dmean$main  = file.path(o.form$main,"qq_dmean")
   if (is.figure && ! file.exists(o.qq.dmean$main)) dir.create(o.qq.dmean$main)
   o.form$qq.dmean  = o.qq.dmean
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the density function of daily means.                             #
   #---------------------------------------------------------------------------------------#
   o.density.dmean       = list()
   o.density.dmean$main  = file.path(o.form$main,"density_dmean")
   if (is.figure && ! file.exists(o.density.dmean$main)) dir.create(o.density.dmean$main)
   o.form$density.dmean  = o.density.dmean
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the skill diagrams.                                              #
   #---------------------------------------------------------------------------------------#
   o.skill                = list()
   o.skill$main           = file.path(o.form$main,"skill")
   if (is.figure && ! file.exists(o.skill$main)) dir.create(o.skill$main)
   for (d in sequence(ndiel)){
      this.diel        = diel.key [d]
      o.diel           = file.path(o.skill$main,this.diel  )
      if (is.figure){
         if (! file.exists(o.diel)) dir.create(o.diel)
      }#end if (is.figure)
      #------------------------------------------------------------------------------------#
      o.skill[[this.diel]] = o.diel
   }#end for (d in 1:ndiel)
   o.form$skill = o.skill
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the Taylor diagrams.                                             #
   #---------------------------------------------------------------------------------------#
   o.taylor      = list()
   o.taylor$main = file.path(o.form$main,"taylor")
   if (is.figure && ! file.exists(o.taylor$main)) dir.create(o.taylor$main)
   for (d in sequence(ndiel)){
      this.diel  = diel.key [d]
      o.diel     = file.path(o.taylor$main,this.diel)
      if (is.figure){
         if (! file.exists(o.diel)) dir.create(o.diel)
      }#end if (is.figure)
      #------------------------------------------------------------------------------------#
      o.taylor[[this.diel]] = o.diel
   }#end for (d in 1:ndiel)
   o.form$taylor = o.taylor
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the light response curve.                                        #
   #---------------------------------------------------------------------------------------#
   o.light       = list()
   o.light$main  = file.path(o.form$main,"light")
   if (is.figure && ! file.exists(o.light$main)) dir.create(o.light$main)
   o.form$light  = o.light
   #---------------------------------------------------------------------------------------#


   #----- Save the full list to the main path list. ---------------------------------------#
   out[[this.form]] = o.form
   #---------------------------------------------------------------------------------------#
}#end for (o in 1:nout)
#------------------------------------------------------------------------------------------#







#------------------------------------------------------------------------------------------#
#      Retrieve all data.                                                                  #
#------------------------------------------------------------------------------------------#
if (reload.hour && file.exists(rdata.hour)){
   cat (" + Reloading hourly data from ",basename(rdata.hour),"...","\n")
   dummy = load(rdata.hour)
}else{
   cat (" + Generating hourly data...","\n")


   #---------------------------------------------------------------------------------------#
   #   Loop through the sites.                                                             #
   #---------------------------------------------------------------------------------------#
   cat (" + Add season and diel keys to seasons and diel...","\n")
   for (p in 1:nsites){
      #----- Grab the observation. --------------------------------------------------------#
      obser      = get(paste("obs",sites$iata[p],sep="."))
      #------------------------------------------------------------------------------------#


      #----- Create some variables to describe season and time of the day. ----------------#
      obser$season    = season(obser$when,add.year=FALSE)
      obser$diel      = 1 + (! obser$nighttime) + obser$highsun
      obser$fortnight = numfortnight  (obser$when)
      obser$toftnight = fnyear.2.chron(obser$when)
      obser$hr.idx    = 3 * floor(obser$hour / 3)
      #------------------------------------------------------------------------------------#


      #----- Fix the measured flags for GPP and ecosystem respiration. --------------------#
      if (all(! obser$measured.cflxst)){
         measured.eco = obser$measured.cflxca
      }else{
         measured.eco = obser$measured.cflxca & obser$measured.cflxst
      }#end if
      obser$measured.nep  = measured.eco
      obser$measured.nee  = measured.eco
      obser$measured.gpp  = measured.eco
      obser$measured.reco = measured.eco & obser$nighttime
      #------------------------------------------------------------------------------------#



      #----- Save the variables to the observations. --------------------------------------#
      dummy = assign(paste("obs",sites$iata[p],sep="."),obser)
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#






   res = list()
   for (p in sequence(nsites)){
      #----- Get the basic information. ---------------------------------------------------#
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
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Get observations, and make sure everything matches.                            #
      #------------------------------------------------------------------------------------#
      obser      = get(paste("obs",iata,sep="."))
      nobser     = length(obser$when)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Get all the statistics and actual values for every simulation.                 #
      #------------------------------------------------------------------------------------#
      trimmed    = FALSE
      for (s in sequence(nsimul)){
         cat("    * Simulation: ",simul$desc[s],"...","\n")

         #----- Load hourly averages. -----------------------------------------------------#
         ans.name = paste("t",iata,"_",simul$name[s],sep="")
         ans.path = file.path(here,ans.name)
         ans.file = file.path(ans.path,"rdata_hour",paste(ans.name,".RData",sep=""))
         load(ans.file)
         nmodel   = length(model$when)
         #---------------------------------------------------------------------------------#


         #=================================================================================#
         #=================================================================================#
         #      In case the models and observations don't match perfectly, we trim         #
         # observations.  This can be done only once, otherwise the models don't have the  #
         # same length, which is unacceptable.                                             #
         #---------------------------------------------------------------------------------#
         if (nobser != nmodel && ! trimmed){
            trimmed = TRUE

            #----- Get the model time range. ----------------------------------------------#
            this.whena = min(as.numeric(model$when))
            this.whenz = max(as.numeric(model$when))
            sel.obser  = obser$when >= this.whena & obser$when <= this.whenz
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all observed variables and trim them to the same model period. #
            #------------------------------------------------------------------------------#
            for ( vname in names(obser)){
               if       (length(obser[[vname]]) == nobser){
                  obser[[vname]] = obser[[vname]][sel.obser]
               }else if (length(obser[[vname]]) != 1     ){
                  sel.now = obser[[vname]] >= this.whena & obser[[vname]] <= this.whenz
                  obser[[vname]] = obser[[vname]][sel.now]
               }#end if
            }#end for
            nobser     = length(obser$when)


            #----- Save the variables to the observations. --------------------------------#
            dummy = assign(paste("obs",sites$iata[p],sep="."),obser)
            #------------------------------------------------------------------------------#
         }else if (nobser == nmodel){
            #----- Datasets match, switch trimmed to TRUE. --------------------------------#
            trimmed = TRUE
            #------------------------------------------------------------------------------#
         }else{
            cat(" -> Simulation:"   ,ans.name          ,"\n")
            cat(" -> Length(obser):",length(obser$when),"\n")
            cat(" -> Length(model):",length(model$when),"\n")
            stop(" Model and obser must have the same length")
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Check whether the model and observations are synchronised.                 #
         #---------------------------------------------------------------------------------#
         if (any(as.numeric(model$when-obser$when) > 1/48,na.rm=TRUE)){
            stop(" All times in the model and observations must match!!!")
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all variables and save daily averages.                           #
         #---------------------------------------------------------------------------------#
         if (s == 1){
            #----- Re-create dates and seasons for daily means. ---------------------------#
            obser$today = dates(obser$when)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Find the daily mean and create flags on whether to use it or not.       #
            #------------------------------------------------------------------------------#
            for (v in sequence(ncompvar)){
               this.compvar  = compvar[[v]]
               this.vnam     = this.compvar$vnam
               this.dmean    = paste("dmean"   ,this.vnam,sep=".")
               this.measured = paste("measured",this.vnam,sep=".")
               this.desc     = this.compvar$desc
               this.unit     = this.compvar$unit
               this.sunvar   = this.compvar$sunvar
               cat("     * ",this.desc,"...","\n")

               #----- Discard all gap-filled entries (except if this is GPP). -------------#
               if (this.vnam == "reco"){
                  discard = ! ( is.finite(obser[[this.vnam]]) & obser$measured.nee )
               }else{
                  discard = ! ( is.finite(obser[[this.vnam]]) & obser[[this.measured]])
               }#end if
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Replace gap-filled values by NA, and find the daily mean.  We replace #
               # the original values by these, so we can eliminate model values easily,    #
               # just adding 0*observations.                                               #
               #---------------------------------------------------------------------------#
               obs.now     = ifelse(discard,NA,obser[[this.vnam]])
               dmean.obser = tapply( X     = obs.now
                                   , INDEX = obser$today
                                   , FUN   = mean
                                   , na.rm = TRUE
                                   )#end tapply
               dmean.count = tapply( X     = is.finite(obs.now)
                                   , INDEX = obser$today
                                   , FUN   = sum
                                   , na.rm = TRUE
                                   )#end tapply
               dmean.obser = ifelse(dmean.count >= nhour.min,dmean.obser,NA)
               #---------------------------------------------------------------------------#



               #----- Save the data. ------------------------------------------------------#
               obser[[this.vnam ]] = obs.now
               obser[[this.dmean]] = dmean.obser
               #---------------------------------------------------------------------------#
            }#end for (v in sequence(ncompvar))
            #------------------------------------------------------------------------------#



            #----- Save the variables to the observations. --------------------------------#
            dummy = assign(paste("obs",sites$iata[p],sep="."),obser)
            #------------------------------------------------------------------------------#
         }#end if (s == 1)
         #---------------------------------------------------------------------------------#





         #----- Create some variables to describe season and time of the day. -------------#
         model$today     = obser$today
         model$season    = obser$season
         model$diel      = obser$diel
         model$fortnight = obser$fortnight
         model$toftnight = obser$toftnight
         model$hr.idx    = obser$hr.idx
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Create daily averages.                                                     #
         #---------------------------------------------------------------------------------#
         dist.comp = list()
         for (v in sequence(ncompvar)){
            this.compvar  = compvar[[v]]
            this.vnam     = this.compvar$vnam
            this.dmean    = paste("dmean"   ,this.vnam,sep=".")
            this.measured = paste("measured",this.vnam,sep=".")
            this.desc     = this.compvar$desc
            this.unit     = this.compvar$unit
            this.sunvar   = this.compvar$sunvar
            cat("     * ",this.desc,"...","\n")

            #----- Discard all gap-filled entries (except if this is GPP). ----------------#
            mod.now = model[[this.vnam]] + 0. * obser[[this.vnam]]
            res.now = obs.now - mod.now
            #------------------------------------------------------------------------------#


            #------ Daily mean of observations. -------------------------------------------#
            dmean.model  = ( tapply(X=mod.now,INDEX=model$today,FUN=mean,na.rm=TRUE)
                           + 0. * obser[[this.dmean]] )
            dmean.today  = chron(names(dmean.model))
            dmean.season = season(dmean.today,add.year=FALSE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Save the variables for later.                                            #
            #------------------------------------------------------------------------------#
            model[[this.vnam ]] = mod.now
            model[[this.dmean]] = dmean.model
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Initialise the list to keep comparisons between model and observations.  #
            #------------------------------------------------------------------------------#
            mat  = matrix( data     = NA
                         , nrow     = ndiel
                         , ncol     = nseason
                         , dimnames = list(diel.key,season.key)
                         )#end matrix
            arr  = array ( data     = NA
                         , dim      = c(ndiel,nseason,nmoment)
                         , dimnames = list(diel.key,season.key,moment.key)
                         )#end array
            comp = list( residuals   = res.now
                       , n           = mat
                       , obs.moment  = arr
                       , mod.moment  = arr
                       , res.moment  = arr
                       , bias        = mat
                       , sigma       = mat
                       , lsq.lnlike  = mat
                       , mse         = mat
                       , rmse        = mat
                       , r.squared   = mat
                       , fvue        = mat
                       , sw.stat     = mat
                       , sw.p.value  = mat
                       , ks.stat     = mat
                       , ks.p.value  = mat
                       )#end list
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Season block.                                                            #
            #------------------------------------------------------------------------------#
            for (ee in sequence(nseason)){
               #---------------------------------------------------------------------------#
               #    Diel block.                                                            #
               #---------------------------------------------------------------------------#
               e.sel = obser$season == ee | ee == nseason
               for (dd in sequence(ndiel)){
                  #------------------------------------------------------------------------#
                  #      Select data.                                                      #
                  #------------------------------------------------------------------------#
                  if (dd == ndiel){
                     #----- Daily means, no daytime check. --------------------------------#
                     sel     = dmean.season == ee | ee == nseason
                     obs.use = obser[[this.dmean]][sel]
                     mod.use = model[[this.dmean]][sel]
                     #---------------------------------------------------------------------#
                  }else{
                     #----- Skip check for sun variables during the night. ----------------#
                     d.sel   = obser$diel == dd | dd == (ndiel-1)
                     s.sel   = obser$highsun | ( ! this.sunvar )
                     sel     = e.sel & d.sel & s.sel
                     obs.use = obser[[this.vnam]][sel]
                     mod.use = model[[this.vnam]][sel]
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #      Skip statistics if no data has been selected or if everything is  #
                  # the same.                                                              #
                  #------------------------------------------------------------------------#
                  sd.obser = sd(obs.use,na.rm=TRUE)
                  #------------------------------------------------------------------------#

                  if (is.finite(sd.obser) && sd.obser > 1.0e-6){
                     #---------------------------------------------------------------------#
                     #      Find the summary of how good (or bad...) the model is.         #
                     #---------------------------------------------------------------------#
                     goodness = test.goodness(x.mod=mod.use,x.obs=obs.use)

                     comp$n         [dd,ee ] = goodness$n
                     comp$obs.moment[dd,ee,] = goodness$obs.moment
                     comp$mod.moment[dd,ee,] = goodness$mod.moment
                     comp$res.moment[dd,ee,] = goodness$res.moment
                     comp$bias      [dd,ee ] = goodness$bias
                     comp$sigma     [dd,ee ] = goodness$sigma
                     comp$lsq.lnlike[dd,ee ] = goodness$lsq.lnlike
                     comp$mse       [dd,ee ] = goodness$mse
                     comp$rmse      [dd,ee ] = goodness$rmse
                     comp$r.squared [dd,ee ] = goodness$r.squared
                     comp$fvue      [dd,ee ] = goodness$fvue
                     comp$sw.stat   [dd,ee ] = goodness$sw.statistic
                     comp$sw.p.value[dd,ee ] = goodness$sw.p.value
                     comp$ks.stat   [dd,ee ] = goodness$ks.statistic
                     comp$ks.p.value[dd,ee ] = goodness$ks.p.value
                     #---------------------------------------------------------------------#
                  }#end if (is.finite(sd.obser) && sd.obser > 1.0e-6)
                  #------------------------------------------------------------------------#
               }#end for (dd in sequence(ndiel))
               #---------------------------------------------------------------------------#
            }#end for (ee in sequence(nseason))
            #------------------------------------------------------------------------------#


            #----- Save the comparison list for this variable. ----------------------------#
            dist.comp[[this.vnam]] = comp
            #------------------------------------------------------------------------------#
         }#end for (v in sequence(ncompvar))
         #---------------------------------------------------------------------------------#



         #----- Save the data and free some memory. ---------------------------------------#
         this$sim[[simul$name[s]]] = dist.comp
         this$ans[[simul$name[s]]] = model
         rm(list=c("dist.comp","model","eddy.complete","eddy.tresume"))
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #----- Copy the data to the results. ------------------------------------------------#
      res[[iata]] = this
      rm(this)
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Save hourly data to RData.                                                       #
   #---------------------------------------------------------------------------------------#
   cat(" + Saving hourly data to ",basename(rdata.hour),"...","\n")
   dummy = save(list=c("res",paste("obs",sites$iata,sep=".")), file=rdata.hour)
   #---------------------------------------------------------------------------------------#
}#end if (reload.hour && file.exists(rdata.hour))
#------------------------------------------------------------------------------------------#















#------------------------------------------------------------------------------------------#
#      Plot the light response curve by site and season.  Each plot has one panel for each #
# simulation.                                                                              #
#------------------------------------------------------------------------------------------#
if (plot.light){
   cat(" + Find the light response curve...","\n")
   for (p in sequence(nsites)){
      iata  = sites$iata[p]
      #----- Get the basic information. ---------------------------------------------------#
      iata          = sites$iata[p]
      this.longname = sites$desc[p]
      cat("   - Site :",this.longname,"...","\n")
      #------------------------------------------------------------------------------------#


      obser = get(paste("obs",iata,sep="."))

      #------------------------------------------------------------------------------------#
      #     Find out when this output variable is finite and measured.  Because PAR and    #
      # total shortwave are very correlated, we ensure that at least one of them is        #
      # measured.                                                                          #
      #------------------------------------------------------------------------------------#
      r.sel = ( ( is.finite(obser$rshort) & obser$measured.rshort )
              | ( is.finite(obser$par)    & obser$measured.par    ) )
      g.sel =   ( is.finite(obser$gpp)    & obser$measured.nee    )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop over the seasons.                                                         #
      #------------------------------------------------------------------------------------#
      for (ee in sequence(nseason)){
         e.sel = obser$season == ee | ee == nseason
         sel   = r.sel & g.sel & e.sel
         n.sel = sum(sel)

         cat("     * ",season.full[ee]," (N = ",n.sel,")...","\n",sep="")


         #---------------------------------------------------------------------------------#
         #     Fit a curve only if there are enough points.                                #
         #---------------------------------------------------------------------------------#
         if (n.sel > 80){
            #------ Observations. ---------------------------------------------------------#
            obs.use = data.frame( par = obser$par[sel]
                                , gpp = obser$gpp[sel] * kgCyr.2.umols
                                )#end data.frame
            #------------------------------------------------------------------------------#



            #----- Select and sort the data. ----------------------------------------------#
            if (tolower(light.method) == "nls"){
               obs.pred = nls.light.response( par.in  = obs.use$par
                                            , gpp     = obs.use$gpp
                                            , first   = c(a1=1,a2=40,a3=500)
                                            , n.boot  = light.n.boot
                                            , control = list(maxiter=500,minFactor=2^-26)
                                            )#end fit.light.response
            }else{
               obs.pred = optim.light.response( par.in  = obs.use$par
                                              , gpp     = obs.use$gpp
                                              , first   = c(a1=1,a2=40,a3=500)
                                              , n.boot  = light.n.boot
                                              , skew    = skew.optim
                                              )#end fit.light.response
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Save range.                                                              #
            #------------------------------------------------------------------------------#
            xlimit = range(obs.use$par,na.rm=TRUE)
            ylimit = range(obs.use$gpp,na.rm=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Grab the data for each simulation.                                       #
            #------------------------------------------------------------------------------#
            mod.use  = list()
            mod.pred = list()
            for (s in sequence(nsimul)){
               mod.now         = res[[iata]]$ans[[s]]
               #----- Select and sort the data. -------------------------------------------#
               mod.use [[s]] = data.frame( par = mod.now$par[sel]
                                         , gpp = mod.now$gpp[sel] * kgCyr.2.umols
                                         )#end data.frame
               if (tolower(light.method) == "nls"){
                  mod.pred[[s]] = nls.light.response( par.in  = mod.use[[s]]$par
                                                    , gpp     = mod.use[[s]]$gpp
                                                    , first   = c(a1=1,a2=40,a3=500)
                                                    , n.boot  = light.n.boot
                                                    , control = list( maxiter   = 500
                                                                    , minFactor = 2^-26
                                                                    )#end list
                                                    )#end fit.light.response
               }else{
                  mod.pred[[s]] = optim.light.response( par.in  = mod.use[[s]]$par
                                                      , gpp     = mod.use[[s]]$gpp
                                                      , first   = c(a1=1,a2=40,a3=500)
                                                      , skew    = skew.optim
                                                      , n.boot  = light.n.boot
                                                      )#end fit.light.response
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Save range.                                                           #
               #---------------------------------------------------------------------------#
               xlimit = range( c(xlimit,mod.use[[s]]$par), na.rm = TRUE )
               ylimit = range( c(ylimit,mod.use[[s]]$gpp), na.rm = TRUE )
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#



            #------ Set some common features. ---------------------------------------------#
            letitre = paste("Light response curve: ",this.longname,"\n"
                           ,season.full[ee]," ( N = ",n.sel,")",sep="")
            lex     = desc.unit(desc="Incoming PAR",unit=untab$umolom2os)
            ley     = desc.unit(desc="Gross Primary Productivity",unit=untab$umolom2os)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Plot the light response curves.                                         #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.light = out[[outform[o]]]$light
               fichier = file.path(out.light,paste("light-",season.suffix[ee],"-",iata,"."
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
                  points( x    = obs.use$par
                        , y    = obs.use$gpp
                        , type = "p"
                        , pch  = 16
                        , col  = grey.bg
                        , cex  = 0.5
                        )#end points
                  points( x    = mod.use[[s]]$par
                        , y    = mod.use[[s]]$gpp
                        , type = "p"
                        , pch  = 16
                        , col  = simul$colour[s]
                        , cex  = 0.5
                        )#end points
                  #------------------------------------------------------------------------#


                  #----- Plot the lines. --------------------------------------------------#
                  lines(x=obs.pred$par,y=obs.pred$gpp ,col=grey.fg,lwd=2,lty="solid" )
                  lines(x=obs.pred$par,y=obs.pred$q025,col=grey.fg,lwd=2,lty="dashed")
                  lines(x=obs.pred$par,y=obs.pred$q975,col=grey.fg,lwd=2,lty="dashed")
                  lines(x   = mod.pred[[s]]$par
                       ,y   = mod.pred[[s]]$gpp
                       ,col = simul$fgcol[s]
                       ,lwd = 2
                       ,lty = "solid"
                       )#end lines
                  lines(x   = mod.pred[[s]]$par
                       ,y   = mod.pred[[s]]$q025
                       ,col = simul$fgcol[s]
                       ,lwd = 2
                       ,lty = "dashed"
                       )#end lines
                  lines(x   = mod.pred[[s]]$par
                       ,y   = mod.pred[[s]]$q975
                       ,col = simul$fgcol[s]
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
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (ee in sequence(nseason))
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#
}#end if (plot.light)
#------------------------------------------------------------------------------------------#




















#------------------------------------------------------------------------------------------#
#      Plot the time series of fortnightly means.                                          #
#------------------------------------------------------------------------------------------#
if (plot.ts.ftnight){
   cat(" + Find the fortnightly means...","\n")
   for (p in sequence(nsites)){
      iata  = sites$iata[p]
      #----- Get the basic information. ---------------------------------------------------#
      iata          = sites$iata[p]
      this.longname = sites$desc[p]
      cat("   - Site :",this.longname,"...","\n")
      #------------------------------------------------------------------------------------#


      obser = get(paste("obs",iata,sep="."))
      today = dates(obser$when)


      #------------------------------------------------------------------------------------#
      #     Loop over variables.                                                           #
      #------------------------------------------------------------------------------------#
      for (v in sequence(ncompvar)){
         this.compvar  = compvar[[v]]
         this.vnam     = this.compvar$vnam
         this.dmean    = paste("dmean"   ,this.vnam,sep=".")
         this.measured = paste("measured",this.vnam,sep=".")
         this.desc     = this.compvar$desc
         this.unit     = this.compvar$unit
         cat("     * ",this.desc,"...","\n")

         #----- Discard all gap-filled entries (except if this is GPP). -------------------#
         obs.now       = obser[[this.vnam]]
         dmean.obser   = obser[[this.dmean]]
         dmean.ftnight = numfortnight(chron(names(dmean.obser)))
         keep          = is.finite(dmean.obser)
         dmean.obser   = dmean.obser  [keep]
         dmean.ftnight = dmean.ftnight[keep]
         #---------------------------------------------------------------------------------#


         #------ Plot only if the data are finite. ----------------------------------------#
         if (length(dmean.obser) > 0){
            #----- Find mean fortnightly period. ------------------------------------------#
            data.in      = data.frame(x=dmean.obser,fortnight=dmean.ftnight)
            boot.out     = boot ( data      = data.in
                                , statistic = boot.fortnight.mean
                                , R         = ftnight.n.boot
                                )#end boot
            obs.expected = apply( X      = boot.out$t
                                , MARGIN = 2
                                , FUN    = mean
                                , na.rm  = TRUE
                                )#end apply
            obs.q025     = apply( X      = boot.out$t
                                , MARGIN = 2
                                , FUN    = quantile
                                , probs  = 0.025
                                , na.rm  = TRUE
                                )#end apply
            obs.q975     = apply( X      = boot.out$t
                                , MARGIN = 2
                                , FUN    = quantile
                                , probs  = 0.975
                                , na.rm  = TRUE
                                )#end apply
            #------------------------------------------------------------------------------#


            #----- Save range. ------------------------------------------------------------#
            ylimit = range(c(obs.expected,obs.q025,obs.q975),na.rm=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            mod.expected = list()
            mod.q025     = list()
            mod.q975     = list()
            for (s in sequence(nsimul)){
               model   = res[[iata]]$ans[[s]]

               #----- Select and sort the data. -------------------------------------------#
               mod.now     = model[[this.vnam ]]
               dmean.model = model[[this.dmean]][keep]
               #---------------------------------------------------------------------------#



               #----- Find mean fortnightly period. ---------------------------------------#
               data.in           = data.frame(x=dmean.model,fortnight=dmean.ftnight)
               boot.out          = boot( data      = data.in
                                       , statistic = boot.fortnight.mean
                                       , R         = ftnight.n.boot
                                       )#end boot
               mod.expected[[s]] = apply( X      = boot.out$t
                                        , MARGIN = 2
                                        , FUN    = mean
                                        , na.rm  = TRUE
                                        )#end apply
               mod.q025    [[s]] = apply( X      = boot.out$t
                                        , MARGIN = 2
                                        , FUN    = quantile
                                        , probs  = 0.025
                                        , na.rm  = TRUE
                                        )#end apply
               mod.q975    [[s]] = apply( X      = boot.out$t
                                        , MARGIN = 2
                                        , FUN    = quantile
                                        , probs  = 0.975
                                        , na.rm  = TRUE
                                        )#end apply
               #---------------------------------------------------------------------------#


               #----- Update range. -------------------------------------------------------#
               ylimit = range(c(ylimit,mod.expected[[s]],mod.q025[[s]],mod.q975[[s]])
                             ,na.rm=TRUE)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#



            #------ Set some common features. ---------------------------------------------#
            letitre = this.longname
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------ Set some common features. ---------------------------------------------#
            xlimit    = c(1,13)
            xat       = c(1:13)
            xlabels   = c(month.abb,month.abb[1])
            x.ftnight = seq(from=1.25,to=12.75,by=0.50)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot the fortnightly means.                                             #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$ts.ftnight[[this.vnam]]
               fichier = file.path(out.now,paste("fnmean-",this.vnam,"-",iata,"."
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
                  #------ Open sub-plot. --------------------------------------------------#
                  par(mar=lo.simul$mar[s,])
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit)
                  if (lo.simul$bottom[s]) axis(side=1,at=xat,labels=xlabels)
                  if (lo.simul$left  [s]) axis(side=2,las=1)
                  abline(v=xat,h=axTicks(2),col=grid.colour,lty="dotted")
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #      Plot the confidence bands (make sure that NA periods are          #
                  # properly skipped).                                                     #
                  #------------------------------------------------------------------------#
                  #----- Find the limits. -------------------------------------------------#
                  x025      = x.ftnight
                  obs.y025  = obs.q025
                  obs.y975  = rev(obs.q975)
                  mod.y025  = mod.q025[[s]]
                  mod.y975  = rev(mod.q975[[s]])
                  keep      = ! (is.na(obs.y025) | is.na(obs.y975))
                  iblck     = cumsum(!keep)
                  #----- Split polygons into finite blocks. -------------------------------#
                  x025      = split (x=x025     [keep],f=iblck[keep])
                  x975      = lapply(X=x025           ,FUN=rev)
                  obs.y025  = split (x=obs.y025 [keep],f=iblck[keep])
                  obs.y975  = split (x=obs.y975 [keep],f=iblck[keep])
                  mod.y025  = split (x=mod.y025 [keep],f=iblck[keep])
                  mod.y975  = split (x=mod.y975 [keep],f=iblck[keep])
                  #----- Plot polygons. ---------------------------------------------------#
                  npoly     = length(x025)
                  for (y in sequence(npoly)){
                     epolygon( x       = c(    x025[[y]],    x975[[y]])
                             , y       = c(obs.y025[[y]],obs.y975[[y]])
                             , col     = grey.bg
                             , angle   = 45
                             , density = 40
                             )#end epolygon
                     epolygon( x       = c(    x025[[y]],    x975[[y]])
                             , y       = c(mod.y025[[y]],mod.y975[[y]])
                             , col     = simul$colour[s]
                             , angle   = -45
                             , density =  40
                             )#end epolygon
                  }#end poly
                  #------------------------------------------------------------------------#



                  #------ Plot the expected values. ---------------------------------------#
                  points( x    = x.ftnight
                        , y    = obs.expected
                        , type = "o"
                        , pch  = 16
                        , col  = grey.fg
                        )#end points
                  points( x    = x.ftnight
                        , y    = mod.expected[[s]]
                        , type = "o"
                        , pch  = 16
                        , col  = simul$fgcol[s]
                        )#end points
                  #------------------------------------------------------------------------#


                  #------ Final stuff. ----------------------------------------------------#
                  box()
                  title(main=simul$desc[s])
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
         }#end if (length(dmean.obser) > 0)
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ncompvar))
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#
}#end if (plot.ts.ftnight)
#------------------------------------------------------------------------------------------#




















#------------------------------------------------------------------------------------------#
#      Plot the time series of fortnightly means.                                          #
#------------------------------------------------------------------------------------------#
if (plot.bp.diel){
   cat(" + Box plot by hour of the day...","\n")
   for (p in sequence(nsites)){
      iata  = sites$iata[p]
      #----- Get the basic information. ---------------------------------------------------#
      iata          = sites$iata[p]
      this.longname = sites$desc[p]
      cat("   - Site :",this.longname,"...","\n")
      #------------------------------------------------------------------------------------#


      obser = get(paste("obs",iata,sep="."))


      #------------------------------------------------------------------------------------#
      #     Loop over variables.                                                           #
      #------------------------------------------------------------------------------------#
      for (v in sequence(ncompvar)){
         this.compvar  = compvar[[v]]
         this.vnam     = this.compvar$vnam
         this.measured = paste("measured",this.vnam,sep=".")
         this.desc     = this.compvar$desc
         this.unit     = this.compvar$unit
         cat("     * ",this.desc,"...","\n")

         #----- Discard all gap-filled entries (except if this is GPP). -------------------#
         discard = ! ( is.finite(obser[[this.vnam]]) & obser[[this.measured]] )
         obs.now = ifelse(discard,NA,obser[[this.vnam]])
         #---------------------------------------------------------------------------------#


         #------ Plot only if the data are finite. ----------------------------------------#
         if (any(is.finite(obs.now))){

            #------ Split data by lists. --------------------------------------------------#
            obs.list        = split(x=obs.now,f=obser$hr.idx)
            names(obs.list) = paste("obs",hour.key,sep=".")
            #------------------------------------------------------------------------------#


            #------ Initialise limits for y. ----------------------------------------------#
            ylimit = range(obs.now,finite=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            mod.list     = list()

            for (s in sequence(nsimul)){
               model   = res[[iata]]$ans[[s]]

               #----- Select and sort the data. -------------------------------------------#
               mod.now = ifelse(discard,NA,model[[this.vnam]])
               #---------------------------------------------------------------------------#


               #----- Split data into lists. ----------------------------------------------#
               mod.list[[s]]        = split(x=mod.now,f=obser$hr.idx)
               names(mod.list[[s]]) = paste("mod",hour.key,sep=".")
               #---------------------------------------------------------------------------#


               #----- Update range. -------------------------------------------------------#
               ylimit  = range(c(ylimit,mod.now),finite=TRUE)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            names(mod.list) = simul$name
            #------------------------------------------------------------------------------#



            #------ Set some common features. ---------------------------------------------#
            letitre = this.longname
            lex     = desc.unit(desc="Time",unit=untab$utc)
            ley     = desc.unit(desc=this.desc,unit=this.unit)
            #------------------------------------------------------------------------------#



            #------ Set some plot defaults. -----------------------------------------------#
            xlimit   = c(0,2*nhour)
            xat      = seq(from=1.5,to=2*nhour-0.5,by=2)
            xlabels  = hour.label
            xgrid    = seq(from=0.5,to=2*nhour+0.5,by=2)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Plot the fortnightly means.                                             #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               out.now = out[[outform[o]]]$bp.diel[[this.vnam]]
               fichier = file.path(out.now,paste("bpdiel-",this.vnam,"-",iata,"."
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
                  bp.list        = c(rbind(obs.list,mod.list[[s]]))
                  names(bp.list) = c(rbind(names(obs.list),names(mod.list[[s]])))
                  bp.colour      = rep(c(grey.bg,simul$col[s]),times=nhour)
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
                  title(main=simul$desc[s])
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
         }#end if (any(is.finite(obs.now)))
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ncompvar))
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
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
      cat("   - ",this.desc,"...","\n")

      #------ Initialise the global Q-Q plot list. ----------------------------------------#
      qq.list = list()
      xlimit  = NULL
      ylimit  = NULL
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Loop over sites.                                                              #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         #----- Grab current site code. ---------------------------------------------------#
         iata            = sites$iata[p]
         #---------------------------------------------------------------------------------#


         #----- Grab data. ----------------------------------------------------------------#
         obser       = get(paste("obs",iata,sep="."))
         dmean.obser = obser[[this.dmean]]
         count.obser = sum(is.finite(dmean.obser))
         #---------------------------------------------------------------------------------#


         #----- Initialise list for this site. --------------------------------------------#
         qq.list[[iata]] = list(n = count.obser)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over simulations.                                                     #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nsimul)){

            #----- Select and sort the data. ----------------------------------------------#
            model       = res[[iata]]$ans[[s]]
            dmean.model = model[[this.dmean]]
            #------------------------------------------------------------------------------#



            #----- Find the Q-Q plot and copy to the global list. -------------------------#
            qq.now = qqplot(x=dmean.obser,y=dmean.model,plot.it=FALSE)
            qq.now = modifyList(x=qq.now,val=list(n=count.obser))
            qq.list[[iata]][[simul$name[s]]] = qq.now
            #------------------------------------------------------------------------------#



            #----- Update range. ----------------------------------------------------------#
            xlimit = range(c(xlimit,qq.now$x),finite=TRUE)
            ylimit = range(c(ylimit,qq.now$y),finite=TRUE)
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------ Set some common features. ---------------------------------------------------#
      letitre = "Q-Q plot of daily means"
      lex     = desc.unit(desc=paste("Observed:" ,this.desc),unit=this.unit)
      ley     = desc.unit(desc=paste("Simulated:",this.desc),unit=this.unit)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the QQ-plots.                                                            #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         out.now = out[[outform[o]]]$qq.dmean
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
            plot.window(xlim=xlimit,ylim=ylimit)
            if (lo.site$bottom[p]) axis.rt(side=1,las=5,off=0.06)
            if (lo.site$left  [p]) axis   (side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            abline(a=0,b=1,col=red.mg,lwd=2,lty="dotdash")
            #------------------------------------------------------------------------------#


            #----- Get Q-Q plot info for this site. ---------------------------------------#
            iata    = sites$iata[p]
            qq.iata = qq.list[[iata]]
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               #------ Plot the lines and points. -----------------------------------------#
               points( x    = qq.iata[[simul$name[s]]]
                     , type = "o"
                     , pch  = 16
                     , cex  = 0.5
                     , col  = simul$col[s]
                     )#end points
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(this.longname,"(N=",qq.iata$n,")",sep=""))
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
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.qq.dmean)
#------------------------------------------------------------------------------------------#




















#------------------------------------------------------------------------------------------#
#      Plot the Q-Q plots of daily means.                                                  #
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
      cat("   - ",this.desc,"...","\n")

      #------ Save all daily means into a list. -------------------------------------------#
      dmean.list = list()
      qlimit     = NULL
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      First loop: obtain the daily means.                                           #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         #----- Initialise list for this site. --------------------------------------------#
         iata               = sites$iata[p]
         #---------------------------------------------------------------------------------#


         #----- Grab data. ----------------------------------------------------------------#
         obser       = get(paste("obs",iata,sep="."))
         dmean.obser = obser[[this.dmean]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      First loop over simulations, to obtain the daily means.                    #
         #---------------------------------------------------------------------------------#
         dmean.model = list()
         for (s in sequence(nsimul)){

            #----- Select and grab the data. ----------------------------------------------#
            model            = res[[iata]]$ans[[s]]
            dmean.model[[s]] = model[[this.dmean]]
            qlimit      = range(c(qlimit,dmean.model[[s]]),finite=TRUE)
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#


         #----- List of daily means. ------------------------------------------------------#
         dmean.list[[iata]] = list( obser = dmean.obser
                                  , model = dmean.model
                                  , n     = sum(is.finite(dmean.obser))
                                  )#end list
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Standardise the quantiles for the density function.  Always choose a power of  #
      # two for the total length, because it reduces the amount of interpolation.          #
      #------------------------------------------------------------------------------------#
      qa    = qlimit[1]
      qz    = qlimit[2]
      quant = seq(from=qa,to=qz,length.out=n.quant)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Second loop: obtain the density functions.                                    #
      #------------------------------------------------------------------------------------#
      dlimit    = NULL
      dens.list = NULL
      for (p in sequence(nsites)){
         #----- Initialise list for this site. --------------------------------------------#
         iata  = sites$iata[p]
         dmean = dmean.list[[iata]]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Check whether the site has data or not.                                    #
         #---------------------------------------------------------------------------------#
         if (any(is.finite(dmean$obser))){


            #----- Density function of observations. --------------------------------------#
            dens.obser = density(x=dmean$obser,n=n.quant,from=qa,to=qz,na.rm=TRUE)$y
            dlimit     = range(c(dlimit,dens.obser),finite=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      First loop over simulations, to obtain the daily means.                 #
            #------------------------------------------------------------------------------#
            dens.model = list()
            for (s in sequence(nsimul)){
               #----- Find daily mean (but only for those days with full record). ---------#
               dens.model[[s]] = density( x     = dmean$model[[s]]
                                        , n     = n.quant
                                        , from  = qa
                                        , to    = qz
                                        , na.rm = TRUE)$y
               dlimit          = range(c(dlimit,dens.model[[s]]),finite=TRUE)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }else{

            #----- Make vector of NAs for density functions. ------------------------------#
            dens.obser = NA + quant
            for (s in sequence(nsimul)){
               dens.model[[s]] = NA + quant
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#
         }#end if (any(is.finite(dmean$obser)))
         #---------------------------------------------------------------------------------#


         #----- List of daily means. ------------------------------------------------------#
         dens.list[[iata]] = list( obser = dens.obser
                                 , model = dens.model
                                 , n     = dmean$n
                                 )#end list
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
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
         out.now = out[[outform[o]]]$density.dmean
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
            plot.window(xlim=qlimit,ylim=dlimit)
            if (lo.site$bottom[p]) axis(side=1,las=1)
            if (lo.site$left  [p]) axis(side=2,las=1)
            grid(col=grid.colour,lty="dotted")
            #------------------------------------------------------------------------------#


            #----- Get Q-Q plot info for this site. ---------------------------------------#
            iata      = sites$iata[p]
            dens.iata = dens.list[[iata]]
            #------------------------------------------------------------------------------#


            #----- Plot the density function of observations. -----------------------------#
            lines ( x    = quant
                  , y    = dens.iata$obser
                  , type = "l"
                  , lwd  = 2.5
                  , col  = foreground
                  )#end points
               #---------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #      Loop over simulations.                                                  #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               #------ Plot the lines. ----------------------------------------------------#
               lines ( x    = quant
                     , y    = dens.iata$model[[s]]
                     , type = "l"
                     , lwd  = 1.6
                     , col  = simul$col[s]
                     )#end points
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------ Final stuff. ----------------------------------------------------------#
            box()
            title(main=paste(this.longname," (N=",dens.iata$n,")",sep=""))
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
         this.measured = paste("measured",this.vnam,sep=".")
         #---------------------------------------------------------------------------------#
         #     Loop over all sites.                                                        #
         #---------------------------------------------------------------------------------#
         for (p in sequence(nsites)){
            iata = sites$iata[p]
            #------------------------------------------------------------------------------#
            #     Grab the data for this simulation.                                       #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               comp.now     = res[[iata]]$sim[[s]][[this.vnam]]
               good.now     = comp.now[[this.good]]
               sfac.now     = 1. + norm.good * ( sqrt(comp.now$obs.moment[,,2]) - 1 )
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
      this.measured = paste("measured",this.vnam,sep=".")
      this.desc     = compvar[[v]]$desc
      this.unit     = compvar[[v]]$unit
      this.sun      = compvar[[v]]$sunvar
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
            iata  = sites$iata[p]
            obs   = get(paste("obs",iata,sep="."))
            nwhen = length(obs$when)

            #------------------------------------------------------------------------------#
            #      Decide which data set to use.                                           #
            #------------------------------------------------------------------------------#
            if (d == ndiel){
               this.obs = obs[[this.dmean]]
               sel      = is.finite(this.obs)
               this.obs = this.obs[sel]
            }else{
               #----- Select this diel (or everything for all day). -----------------------#
               d.sel    = (obs$diel == d | d == (ndiel-1))
               s.sel    = obs$highsun | (! this.sun)
               o.sel    = is.finite(obs[[this.vnam]])
               sel      = d.sel & s.sel & o.sel
               sel      = ifelse(is.na(sel),FALSE,sel)
               this.obs = obs[[this.vnam]][sel]
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Find the standard deviation of this observation.  Skip the site if       #
            # everything is zero.                                                          #
            #------------------------------------------------------------------------------#
            sdev.obs.now = sd(this.obs,na.rm=TRUE)
            sel          = sel & is.finite(sdev.obs.now) & sdev.obs.now > 0
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
                  #      Decide which data set to use.                                     #
                  #------------------------------------------------------------------------#
                  mod  = res[[iata]]$ans[[simul.key[s]]]
                  if (d == ndiel){
                     this.mod = mod[[this.dmean]][sel]
                  }else{
                     #----- Copy simulation. -------------------------------------------------#
                     this.mod = mod[[this.vnam]][sel]
                     #------------------------------------------------------------------------#
                  }#end if
                  mod.diel[[iata]][,s] = this.mod
                  #------------------------------------------------------------------------#


                  #----- Check number of valid entries. -----------------------------------#
                  if (! is.null(cnt.diel[[iata]])){
                     cnt.diel[[iata]] = sum(is.finite(this.obs))
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Find the normalised bias and model standard deviation. -----------#
                  comp        = res[[iata]]$sim[[simul.key[s]]][[this.vnam]]
                  bias.now    = comp$bias [d,nseason] / sqrt(comp$obs.moment[d,nseason,2])
                  sigma.now   = comp$sigma[d,nseason] / sqrt(comp$obs.moment[d,nseason,2])
                  bias.range  = c(bias.range ,bias.now   )
                  sigma.range = c(sigma.range,sigma.now  )
                  #------------------------------------------------------------------------#
               }#end for (s in sequence(nsimul))
               #---------------------------------------------------------------------------#
            }#end if (any(sel))
            #------------------------------------------------------------------------------#
         }#end for (p in 1:nsites)
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
            }#end for (o in 1:nout) 
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         }#end if (length(ref) > 2)
         #---------------------------------------------------------------------------------#
      }#end for (d in 1:ndiel)
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:ncompvar)
   #---------------------------------------------------------------------------------------#
}#end if (plot.skill.taylor)
#------------------------------------------------------------------------------------------#

