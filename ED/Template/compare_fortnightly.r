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
outroot = file.path(here,paste("fortnightly_comp_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#



#----- Info on fortnightly data. ----------------------------------------------------------#
reload.fortnight = TRUE
rdata.path       = file.path(here,"RData_fortnightly")
rdata.fortnight  = file.path(rdata.path,"fortnightly_ed22.RData")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Site settings:                                                                       #
# eort      -- first letter ("e" or "t")                                                   #
# sites     -- site codes ("IATA")                                                         #
# sites.pch -- site symbols                                                                #
#------------------------------------------------------------------------------------------#
eort           = "t"
sites          = c("gyf","s67","s83","pdg","rja","m34")  # ,"pnz","ban","cax"
sites.pch      = c(    2,    5,    9,   13,    1,    6)  # ,    4,    8,    0
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
                     , desc     = c("Big leaf, 2 PFTs"    ,"Big leaf, 5 PFTs"
                                   ,"Size only, 2 PFTs"   ,"Size only, 5 PFTs"
                                   ,"Size and age, 2 PFTs","Size and age, 5 PFTs"
                                   )#end c
                     , verbose  = c("Big leaf, 2 PFTs"    ,"Big leaf, 5 PFTs"
                                   ,"Size only, 2 PFTs"   ,"Size only, 5 PFTs"
                                   ,"Size and age, 2 PFTs","Size and age, 5 PFTs"
                                   )#end c
                     , colour   = c("slateblue4" ,"purple1"     
                                   ,"dodgerblue3","deepskyblue" 
                                   ,"chartreuse4","chartreuse"  
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
ptsz           = 18                    # Font size.
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

min.fortnight  = 1/3                   # Minimum number of observations to consider the
                                       #     fortnight period.
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
simlty.key     = rep("solid",times=n.sim)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Dump the information to a list.                                                      #
#------------------------------------------------------------------------------------------#
simul       = data.frame( name             = simul.key
                        , desc             = simleg.key
                        , colour           = simcol.key
                        , lty              = simlty.key
                        , stringsAsFactors = FALSE
                        )#end data.frame
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      List the keys for all dimensions.                                                   #
#------------------------------------------------------------------------------------------#
sites.key   = sites
sites.desc  = poilist$longname[match(sites.key,poilist$iata)]
control.key = apply(X = sapply(X=control,FUN=c),MARGIN=1,FUN=unlist)[,"vnam"]
compvar.key = apply(X = sapply(X=compvar,FUN=c),MARGIN=1,FUN=unlist)$vnam
compvar.sym = apply(X = sapply(X=compvar,FUN=c),MARGIN=1,FUN=unlist)$symbol
good.key    = apply(X = sapply(X=good   ,FUN=c),MARGIN=1,FUN=unlist)[,"vnam"]

hour.num    = seq(from=0,to=23,by=3)
diel.key    = c( paste( sprintf("%2.2i",(hour.num - 1) %% 24)
                      , sprintf("%2.2i",(hour.num + 2) %% 24)
                      , sep = "."
                      )#end paste
               , "mean.night", "mean.twilight", "mean.day", "fortnight.mean"
               , "all.night" , "all.twilight" , "all.day" , "all.hrs"
               )#end c
diel.desc   = c( paste( sprintf("%2.2i",(hour.num - 1) %% 24)
                      , "-"
                      , sprintf("%2.2i",(hour.num + 2) %% 24)
                      , " UTC"
                      , sep = ""
                      )#end paste
               , "Nocturnal mean", "Twilight mean", "Diurnal mean", "Fortnightly mean"
               , "Nighttime"     , "Twilight"     , "Daytime"     , "All hours"
               )#end c
#------------------------------------------------------------------------------------------#



#----- Set the various dimensions associated with variables, simulations, and sites. ------#
nsites   = length(sites.key  )
nsimul   = length(simul.key  )
ncompvar = length(compvar.key)
ncontrol = length(control.key)
ngood    = length(good.key   )
ndiel    = length(diel.key   )
#------------------------------------------------------------------------------------------#


#----- Set the time periods for which we aggregate the data. ------------------------------#
naggr     = ndiel-4
aggr.key  = diel.key [sequence(naggr)]
aggr.desc = diel.desc[sequence(naggr)]
#------------------------------------------------------------------------------------------#




#----- Load observations. -----------------------------------------------------------------#
obser.file = paste(srcdir,"LBA_MIP.nogapfill.RData",sep="/")
load(file=obser.file)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
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


   #----- Save the full list to the main path list. ---------------------------------------#
   out[[this.form]] = o.form
   #---------------------------------------------------------------------------------------#
}#end for (o in 1:nout)
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#   Loop through the sites.                                                                #
#------------------------------------------------------------------------------------------#
if (reload.fortnight && file.exists(rdata.fortnight)){
   cat (" + Reloading fortnightly data from ",basename(rdata.fortnight),"...","\n")
   dummy = load(rdata.fortnight)
}else{
   cat (" + Generating fortnightly data...","\n")


   cat (" + Add fortnight and diel keys to observations...","\n")
   for (p in sequence(nsites)){
      #----- Grab the observation. --------------------------------------------------------#
      obser = get(paste("obs",sites[p],sep="."))
      #------------------------------------------------------------------------------------#


      #----- Create some variables to describe season and time of the day. ----------------#
      if (! "fortnight" %in% names(obser)){
         obser$fortnight = numfortnight(obser$when)
      }#end if
      if (! "diel"      %in% names(obser)){
         obser$diel      = (! obser$nighttime) + obser$highsun
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Save the variables to the observations. --------------------------------------#
      dummy = assign(paste("obs",sites[p],sep="."),obser)
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Retrieve all data.                                                               #
   #---------------------------------------------------------------------------------------#
   cat (" + Retrieve model results for all sites...","\n")
   res = list()
   for (p in sequence(nsites)){
      #----- Get the basic information. ---------------------------------------------------#
      iata          = sites[p]
      im            = match(iata,poilist$iata)
      this          = list()
      this$short    = poilist$short   [im]
      this$longname = poilist$longname[im]
      this$iata     = poilist$iata    [im]
      this$lon      = poilist$lon     [im]
      this$lat      = poilist$lat     [im]
      this$sim      = list()
      cat("   - Site :",this$longname,"...","\n")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Get all the statistics and actual values for every simulation.                 #
      #------------------------------------------------------------------------------------#
      trimmed = FALSE
      for (s in sequence(nsimul)){
         cat("    * Simulation: ",simul$desc[s],"...","\n")
         ans.name = paste(eort,iata,"_",simul$name[s],sep="")
         ans.path = file.path(here,ans.name)
         ans.file = paste(ans.path,"rdata_hour",paste(ans.name,".RData",sep="")
                         ,sep="/")
         load(ans.file)


         #---------------------------------------------------------------------------------#
         #     Get observations, and make sure everything matches.                         #
         #---------------------------------------------------------------------------------#
         obser      = get(paste("obs",iata,sep="."))

         nobser     = length(obser$when)
         nmodel     = length(model$when)


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
            #------------------------------------------------------------------------------#


            #----- Save the variables to the observations. --------------------------------#
            dummy = assign(paste("obs",sites[p],sep="."),obser)
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



         #------ Create indices for 3-hour blocks and fortnight periods. ------------------#
         obser$toftnight = fnyear.2.chron(obser$when)
         obser$hr.idx    = 3*floor(obser$hour        / 3)
         model$toftnight = fnyear.2.chron(model$when)
         model$hr.idx    = 3*floor(hours(model$when) / 3)
         #---------------------------------------------------------------------------------#




         #------ Define empty templates with readable names to host the data. -------------#
         fortnight.key    = chron(sort(unique(obser$toftnight)))
         nfortnight       = length(fortnight.key)
         empty.aggr       = array( data     = NA
                                 , dim      = c(nfortnight,naggr)
                                 , dimnames = list(paste(fortnight.key),aggr.key)
                                 )#end array
         empty.vec        = rep   ( NA,times=ndiel)
         names(empty.vec) = diel.key
         empty.mat        = matrix( data     = NA
                                  , nrow     = ndiel
                                  , ncol     = 4
                                  , dimnames = list( diel.key
                                                   , c("mean","variance"
                                                      ,"skewness","kurtosis")
                                                   )#end list
                                  )#end matrix
         #---------------------------------------------------------------------------------#



         #----- Template for each variable. -----------------------------------------------#
         var.template  = list( model        = empty.aggr
                             , obser        = empty.aggr
                             , resid        = empty.aggr
                             , count        = empty.aggr
                             , total        = empty.aggr
                             , bias         = empty.vec
                             , sigma        = empty.vec
                             , rmse         = empty.vec
                             , r.squared    = empty.vec
                             , fvue         = empty.vec
                             , sw.stat      = empty.vec
                             , sw.p.value   = empty.vec
                             , ks.stat      = empty.vec
                             , ks.p.value   = empty.vec
                             , lsq.lnlike   = empty.vec
                             , obser.moment = empty.mat
                             , model.moment = empty.mat
                             , resid.moment = empty.mat
                             )#end list
         #---------------------------------------------------------------------------------#



         #----- Initialise list that will hold 
         vsumm         = list( fortnight    = sort(unique(obser$toftnight))
                             , hour         = sort(unique(obser$hour     ))
                             , diel         = empty.aggr
                             , total        = empty.aggr
                             )#end list
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Define the hours status regarding the diel.                                #
         #---------------------------------------------------------------------------------#
         vsumm$diel [,] = cbind( tapply( X     = obser$diel
                                       , INDEX = list(obser$toftnight,obser$hr.idx)
                                       , FUN   = commonest
                                       , na.rm = TRUE
                                       )#end tapply
                               , tapply( X     = obser$diel + NA
                                       , INDEX = list(obser$toftnight,obser$diel)
                                       , FUN   = commonest
                                       , na.rm = TRUE
                                       )#end tapply
                               , tapply( X     = obser$diel + NA
                                       , INDEX = obser$toftnight
                                       , FUN   = commonest
                                       , na.rm = TRUE
                                       )#end tapply
                               )#end cbind
         vsumm$diel[! is.finite(vsumm$diel)] = NA

         vsumm$total[,] = cbind( tapply( X     = is.finite(obser$when)
                                       , INDEX = list(obser$toftnight,obser$hr.idx)
                                       , FUN   = sum
                                       , na.rm = TRUE
                                       )#end tapply
                               , tapply( X     = is.finite(obser$when)
                                       , INDEX = list(obser$toftnight,obser$diel)
                                       , FUN   = sum
                                       , na.rm = TRUE
                                       )#end tapply
                               , tapply( X     = is.finite(obser$when)
                                       , INDEX = obser$toftnight
                                       , FUN   = sum
                                       , na.rm = TRUE
                                       )#end tapply
                               )#end cbind
         vsumm$total[! is.finite(vsumm$total)] = 0
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Loop over all variables, and build the array with information.              #
         #---------------------------------------------------------------------------------#
         for (v in sequence(ncompvar)){
            this.comp     = compvar[[v]]
            this.vnam     = this.comp$vnam
            this.measured = paste("measured",this.vnam,sep=".")
            this.symbol   = this.comp$symbol
            this.desc     = this.comp$desc
            this.unit     = this.comp$unit
            col.obser     = this.comp$col.obser
            col.model     = this.comp$col.model
            leg.corner    = this.comp$leg.corner
            sunvar        = this.comp$sunvar



            #------ Use now as the place holder. ------------------------------------------#
            ans.now     = var.template
            model.now   = model[[this.vnam    ]]
            obser.now   = obser[[this.vnam    ]]
            measr.now   = obser[[this.measured]]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Correct sel in case this is a "Sun" variable.  We only check Sun         #
            # variables during daytime.                                                    #
            #------------------------------------------------------------------------------#
            if (sunvar){
               sel.sun = obser$daytime
            }else{
               sel.sun = rep(TRUE,times=length(obser$daytime))
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Use only the model results from times when measurements existed, so we   #
            # won't bias the model.                                                        #
            #------------------------------------------------------------------------------#
            skip            = ! ( sel.sun & measr.now)
            model.now[skip] = NA
            obser.now[skip] = NA
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Aggregate model and observations.                                       #
            #------------------------------------------------------------------------------#
            ans.now$model[,] = cbind( tapply( X     = model.now
                                            , INDEX = list(obser$toftnight,obser$hr.idx)
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = model.now
                                            , INDEX = list(obser$toftnight,obser$diel)
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = model.now
                                            , INDEX = obser$toftnight
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    )#end cbind
            ans.now$obser[,] = cbind( tapply( X     = obser.now
                                            , INDEX = list(obser$toftnight,obser$hr.idx)
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = obser.now
                                            , INDEX = list(obser$toftnight,obser$diel)
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = obser.now
                                            , INDEX = obser$toftnight
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    )#end cbind

            ans.now$resid[,] = cbind( tapply( X     = model.now - obser.now
                                            , INDEX = list(obser$toftnight,obser$hr.idx)
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = model.now - obser.now
                                            , INDEX = list(obser$toftnight,obser$diel)
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = model.now - obser.now
                                            , INDEX = obser$toftnight
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            )#end tapply
                                    )#end cbind

            ans.now$count[,] = cbind( tapply( X     = is.finite(model.now)
                                            , INDEX = list(obser$toftnight,obser$hr.idx)
                                            , FUN   = sum
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = is.finite(model.now)
                                            , INDEX = list(obser$toftnight,obser$diel)
                                            , FUN   = sum
                                            , na.rm = TRUE
                                            )#end tapply
                                    , tapply( X     = is.finite(model.now)
                                            , INDEX = obser$toftnight
                                            , FUN   = sum
                                            , na.rm = TRUE
                                            )#end tapply
                                    )#end cbind
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #       Replace all "junk" by NA.                                              #
            #------------------------------------------------------------------------------#
            ans.now$model[! is.finite(ans.now$model)] = NA
            ans.now$obser[! is.finite(ans.now$obser)] = NA
            ans.now$resid[! is.finite(ans.now$resid)] = NA
            ans.now$count[! is.finite(ans.now$count)] = 0
            
            sparse                = ans.now$count / vsumm$total < min.fortnight
            sparse                = sparse | vsumm$total == 0
            ans.now$model[sparse] = NA
            ans.now$obser[sparse] = NA
            ans.now$resid[sparse] = NA
            ans.now$count[sparse] = 0
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #       Find the statistics for rows.                                          #
            #------------------------------------------------------------------------------#
            for (dh in sequence(ndiel)){
               #---------------------------------------------------------------------------#
               #     Select the period of the day to find statistics.                      #
               #---------------------------------------------------------------------------#
               if (dh <= naggr){
                  this.obser = ans.now$obser[,dh]
                  this.model = ans.now$model[,dh]
                  this.resid = ans.now$resid[,dh]
               }else{
                  sel        = ( is.finite(vsumm$diel)
                               & (dh == ndiel | vsumm$diel == (dh-naggr-1) )
                               )#end sel
                  this.obser = ans.now$obser[sel]
                  this.model = ans.now$model[sel]
                  this.resid = ans.now$resid[sel]
               }#end if
               #---------------------------------------------------------------------------#



               #----- Find and plot the distribution function for this time. --------------#
               sd.obser = sd(this.obser,na.rm=TRUE)
               if (is.finite(sd.obser) && sd.obser > 1.0e-6){


                  #------------------------------------------------------------------------#
                  #      Find multiple statistics that may be used for finding the support #
                  # function.                                                              #
                  #------------------------------------------------------------------------#
                  o.sn.stats = sn.stats   (this.obser,na.rm=TRUE)
                  o.location = o.sn.stats[1]
                  o.scale    = o.sn.stats[2]
                  o.shape    = o.sn.stats[3]
                  o.mean     = mean       (this.obser,na.rm=TRUE)
                  o.vari     = var        (this.obser,na.rm=TRUE)
                  o.sdev     = sd         (this.obser,na.rm=TRUE)
                  o.skew     = skew       (this.obser,na.rm=TRUE)
                  o.kurt     = kurt       (this.obser,na.rm=TRUE)
                  m.sn.stats = sn.stats   (this.model,na.rm=TRUE)
                  m.location = m.sn.stats[1]
                  m.scale    = m.sn.stats[2]
                  m.shape    = m.sn.stats[3]
                  m.mean     = mean       (this.model,na.rm=TRUE)
                  m.vari     = var        (this.model,na.rm=TRUE)
                  m.sdev     = sd         (this.model,na.rm=TRUE)
                  m.skew     = skew       (this.model,na.rm=TRUE)
                  m.kurt     = kurt       (this.model,na.rm=TRUE)
                  r.sn.stats = sn.stats   (this.resid,na.rm=TRUE)
                  r.location = r.sn.stats[1]
                  r.scale    = r.sn.stats[2]
                  r.shape    = r.sn.stats[3]
                  r.mean     = mean       (this.resid,na.rm=TRUE)
                  r.vari     = var        (this.resid,na.rm=TRUE)
                  r.sdev     = sd         (this.resid,na.rm=TRUE)
                  r.skew     = skew       (this.resid,na.rm=TRUE)
                  r.kurt     = kurt       (this.resid,na.rm=TRUE)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Run a Kolgomorov-Smirnov test comparing the two distributions.    #
                  #------------------------------------------------------------------------#
                  this.ks                = ks.test(x=this.obser,y=this.model)
                  ans.now$ks.stat   [dh] = this.ks$statistic
                  ans.now$ks.p.value[dh] = this.ks$p.value
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Find the mean bias, the root mean square error, the coefficient   #
                  # of determination, and the fraction of variance unexplained for this    #
                  # simulation.                                                            #
                  #------------------------------------------------------------------------#
                  goodness = test.goodness ( x.mod        = this.model
                                           , x.obs        = this.obser
                                           , n.parameters = NULL
                                           )#end test.goodness
                  ans.now$bias      [dh] = goodness$bias
                  ans.now$sigma     [dh] = goodness$sigma
                  ans.now$rmse      [dh] = goodness$rmse
                  ans.now$lsq.lnlike[dh] = goodness$lsq.lnlike
                  ans.now$r.squared [dh] = goodness$r.squared
                  ans.now$fvue      [dh] = goodness$fvue
                  ans.now$sw.stat   [dh] = goodness$sw.statistic
                  ans.now$sw.pvalue [dh] = goodness$sw.pvalue
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #      Find the first four moments of the distribution for observations, #
                  # model, and residuals.                                                  #
                  #------------------------------------------------------------------------#
                  ans.now$obser.moment[dh,] = c(o.mean,o.vari,o.skew,o.kurt)
                  ans.now$model.moment[dh,] = c(m.mean,m.vari,m.skew,m.kurt)
                  ans.now$resid.moment[dh,] = c(r.mean,r.vari,r.skew,r.kurt)
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#
            vsumm[[this.vnam]] = ans.now
         }#end for
         #---------------------------------------------------------------------------------#


         this$sim[[simul$name[s]]] = vsumm
         rm(list=c("obser","model","eddy.complete","eddy.tresume"))
      }#end for
      #------------------------------------------------------------------------------------#


      #----- Copy the data to the results. ------------------------------------------------#
      res[[iata]] = this
      rm(this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Save fortnightly data to RData.                                                  #
   #---------------------------------------------------------------------------------------#
   cat(" + Saving fortnightly data to ",basename(rdata.fortnight),"...","\n")
   dummy = save(res, file=rdata.fortnight)
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Plot the spider web for all sites and all variables, for statistics that may can be #
# plotted in a spider web (i.e., positive defined).                                        #
#------------------------------------------------------------------------------------------#
cat(" + Plot the spider web diagrams for all sites and all variables...","\n")
good.loop = which(unlist(sapply(X=good,FUN=c)["spider",]))
for (g in good.loop){
   #---- Copy structured variables to convenient scratch scalars. -------------------------#
   this.good = good[[g]]$vnam
   desc.good = good[[g]]$desc
   norm.good = good[[g]]$normalise
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Load all data to an array.                                                       #
   #---------------------------------------------------------------------------------------#
   cat("   - ",desc.good,"...","\n")
   web = array( dim      = c   (nsimul,nsites,ncompvar,ndiel)
              , dimnames = list(simul.key,sites.key,compvar.key,diel.key)
              )#end array 
   for (v in sequence(ncompvar)){
      this.vnam     = compvar[[v]]$vnam
      this.measured = paste("measured",this.vnam,sep=".")

      #------------------------------------------------------------------------------------#
      #     Loop over all sites.                                                           #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         iata = sites.key[p]



         #---------------------------------------------------------------------------------#
         #     Grab the data for this simulation.                                          #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nsimul)){
            this         = res[[iata]]$sim[[s]][[this.vnam]]
            good.now     = this[[this.good]]
            if (norm.good){
               sfac = sqrt(this$obser.moment[,2])
               sfac = ifelse(sfac == 0.,NA,1/sfac)
            }else{
               sfac = 1. + 0. * good.now
            }#end if
            web[s,p,v,] = abs(good.now) * sfac
         }#end for (s in 1:nsimul)
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Plot the spider webs by diel and variable.                                        #
   #---------------------------------------------------------------------------------------#
   for (d in sequence(ndiel)){
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Webs by site (all variables).                                                  #
      #------------------------------------------------------------------------------------#
      for (p in sequence(nsites)){
         iata     = sites.key[p]

         letitre = paste(desc.good," - ",sites.desc[p],"\n",diel.desc[d],sep="")

         if (any(is.finite(web[,p,,d]))){
            v.sel = is.finite(colSums(web[,p,,d]))

            if (this.good %in% "sw.stat"){
               web.range = c(0,1)
            }else{
               web.range = range(c(0,web[,p,v.sel,d]),na.rm=TRUE)
            }#end if
            if (ptsz <= 11){
               web.lim   = pretty(web.range,n=5)
            }else if (ptsz <= 14){
               web.lim   = pretty(web.range,n=4)
            }else{
               web.lim   = pretty(web.range,n=3)
            }#end if

            #------------------------------------------------------------------------------#
            #     Webs by variable (all sites).                                            #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into 3, and add site and simulation legends at the   #
               # bottom.                                                                   #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,3.0,0))
               layout(mat = rbind(2,1),height = c(18,4))
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
                      , cex     = 0.8 * cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the spider web.                                                  #
               #---------------------------------------------------------------------------#
               radial.flex( lengths          = web[,p,v.sel,d]
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
         }#end if (any(is.finite(web[,,v,d])))
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:ncompvar)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Webs by variable (all sites).                                                  #
      #------------------------------------------------------------------------------------#
      for (v in sequence(ncompvar)){
         this.vnam     = compvar[[v]]$vnam
         this.desc     = compvar[[v]]$desc

         letitre = paste(desc.good," - ",this.desc,"\n",diel.desc[d],sep="")

         if (any(is.finite(web[,,v,d]))){
            p.sel = is.finite(colSums(web[,,v,d]))

            if (this.good %in% "sw.stat"){
               web.range = c(0,1)
            }else{
               web.range = range(c(0,web[,p.sel,v,d]),na.rm=TRUE)
            }#end if
            web.lim   = pretty(web.range,n=4)

            #------------------------------------------------------------------------------#
            #     Webs by variable (all sites).                                            #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into 3, and add site and simulation legends at the   #
               # bottom.                                                                   #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,3.0,0))
               layout(mat = rbind(2,1),height = c(5.0,1.0))
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
                      , cex     = 0.8 * cex.ptsz
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the spider web.                                                  #
               #---------------------------------------------------------------------------#
               radial.flex( lengths          = web[,p.sel,v,d]
                          , labels           = toupper(sites.key[p.sel])
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
            #------------------------------------------------------------------------------#
         }#end if (any(is.finite(web[,,v,d])))
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ncompvar))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   }#end for (d in sequence(ndiel))
   #---------------------------------------------------------------------------------------#
}#end for (g in good.loop)
#------------------------------------------------------------------------------------------#







#------------------------------------------------------------------------------------------#
#         Plot the various statistics as functions of the site "completion".               #
#------------------------------------------------------------------------------------------#
cat (" + Plot cross-model and cross-site diagrams...","\n")
performance = array( data     = NA
                   , dim      = c(ncompvar,ngood,nsimul)
                   , dimnames = list(compvar.key,good.key,simul.key)
                   )#end array
for (v in sequence(ncompvar)){
   #----- Copy the variable information. --------------------------------------------------#
   this.vnam     = compvar[[v]]$vnam
   this.desc     = compvar[[v]]$desc
   this.unit     = compvar[[v]]$unit
   this.sun      = compvar[[v]]$sunvar
   this.measured = paste("measured",this.vnam,sep=".")
   cat("   - ",this.desc,"...","\n")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Loop over all parts of the day.                                                  #
   #---------------------------------------------------------------------------------------#
   for (dh in sequence(ndiel)){
      cat("     * ",diel.desc[dh],"...","\n")

      #------------------------------------------------------------------------------------#
      #      Loop over all sites, normalise the data and create the vector for the model.  #
      #------------------------------------------------------------------------------------#
      obs.diel    = list()
      mod.diel    = list()
      cnt.diel    = rep(x=NA,times=nsites); names(cnt.diel) = sites
      bias.range  = NULL
      sigma.range = NULL
      for (p in sequence(nsites)){
         iata  = sites[p]
         #----- Load data from this site. -------------------------------------------------#
         this       = res[[iata]]
         #---------------------------------------------------------------------------------#

         #----- Copy the modelled and observed data, and update ranges. -------------------#
         for (s in sequence(nsimul)){
            #----- Get the summary of this simulation and load the variable. --------------#
            vsumm   = this$sim[[simul$name[s]]]
            ans.now = vsumm[[this.vnam]]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make sure to initialise the data set summary.                            #
            #------------------------------------------------------------------------------#
            if (s == 1){
               #---------------------------------------------------------------------------#
               #      Matrix size depends upon the type of aggregation.                    #
               #---------------------------------------------------------------------------#
               if (dh <= naggr){
                  fortnight        = vsumm$fortnight
                  nfortnight       = length(fortnight)
                  obs.diel[[iata]] = ans.now$obser[,dh]
                  cnt.diel[[iata]] = sum(ans.now$count[,dh])
                  mod.diel[[iata]] = matrix( ncol     = nsimul
                                           , nrow     = nfortnight
                                           , dimnames = list(paste(fortnight),simul.key)
                                           )#end matrix
               }else{
                  dft.sel          = ( is.finite(vsumm$diel)
                                     & (dh == ndiel | vsumm$diel == (dh-naggr-1) )
                                     )#end sel
                  ndft.sel         = sum(dft.sel)
                  obs.diel[[iata]] = ans.now$obser[dft.sel]
                  cnt.diel[[iata]] = sum(ans.now$count[dft.sel])
                  mod.diel[[iata]] = matrix( ncol     = nsimul
                                           , nrow     = ndft.sel
                                           , dimnames = list(NULL,simul.key)
                                           )#end matrix
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Select the period of the day to find statistics.                         #
            #------------------------------------------------------------------------------#
            if (dh <= naggr){
               mod.diel[[iata]][,s] = ans.now$model[,dh]
            }else{
               mod.diel[[iata]][,s] = ans.now$model[dft.sel]
            }#end if
            #------------------------------------------------------------------------------#



            #----- Find the normalised bias and model standard deviation. -----------------#
            bias.now    = ans.now$bias [dh] / sqrt(ans.now$obser.moment[dh,2])
            sigma.now   = ans.now$sigma[dh] / sqrt(ans.now$obser.moment[dh,2])
            bias.range  = c(bias.range ,bias.now   )
            sigma.range = c(sigma.range,sigma.now  )
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (p in sequence(nsites))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot Taylor and skill plots only if there is anything to plot.                 #
      #------------------------------------------------------------------------------------#
      ok.taylor.skill = ( length(unlist(obs.diel)) > 0  && any(cnt.diel > 0)           &&
                          any(is.finite(bias.range))    && any(is.finite(sigma.range))    )

      if (ok.taylor.skill){
         #---- Fix ranges. ----------------------------------------------------------------#
         xy.range    = 1.04 * max(abs(c(bias.range,sigma.range)),na.rm=TRUE)
         bias.range  = 1.04 * xy.range  * c(-1,1)
         sigma.range = 1.04 * xy.range  * c( 1,0)
         r2.range    = range(1-xy.range^2,1)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Calculate the size for the points in the Skill and Taylor diagrams.  Make   #
         # it proportional to the number of points used to evaluate each place.            #
         #---------------------------------------------------------------------------------#
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
         #---------------------------------------------------------------------------------#



         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Skill plot.                                                                 #
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Plot title.                                                                 #
         #---------------------------------------------------------------------------------#
         letitre = paste(" Skill diagram - ",this.desc,"\n",diel.desc[dh],sep="")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.skill = out[[outform[o]]]$skill[[diel.key[dh]]]
            fichier   = file.path(out.skill,paste("skill-",this.vnam,"-",diel.key[dh]
                                                 ,".",outform[o],sep="")
                                 )#end file.path
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
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into 4, and add site and simulation legends at the      #
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
                   , legend  = toupper(sites.key)
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = sites.pch
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
            #------------------------------------------------------------------------------#




            #----- Legend: the simulations. -----------------------------------------------#
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
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over sites.                                                         #
            #------------------------------------------------------------------------------#
            myskill = NULL
            for (p in sequence(nsites)){
               iata = sites[p]

               #----- Skip the site if there is no data. ----------------------------------#
               ok.iata = length(obs.diel[[iata]]) > 0 && any(is.finite(obs.diel[[iata]]))
               ok.iata = ok.iata && ( ! is.na(ok.iata))
               if (ok.iata){
                  #------------------------------------------------------------------------#
                  #     We call skill twice for each site in case the site has two PCHs.   #
                  #------------------------------------------------------------------------#
                  myskill = skill.plot( obs           = obs.diel[[iata]]
                                      , obs.options   = list( col = foreground
                                                            , cex = 2.0
                                                            )#end list
                                      , mod           = mod.diel[[iata]]
                                      , mod.options   = list( col = simul$col
                                                            , bg  = simul$col
                                                            , pch = sites.pch[p]
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
                  #------------------------------------------------------------------------#
               }#end if (length(obs.diel[[iata]] > 0)
               #---------------------------------------------------------------------------#
            }#end for (p in 1:nsites)
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
         letitre = paste(" Taylor diagram - ",this.desc,"\n",diel.desc[d],sep="")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.taylor = out[[outform[o]]]$taylor[[diel.key[dh]]]
            fichier    = file.path(out.taylor,paste("taylor-",this.vnam,"-",diel.key[d]
                                                   ,".",outform[o],sep=""))
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
                   , legend  = toupper(sites.key)
                   , col     = foreground
                   , pt.bg   = foreground
                   , pch     = sites.pch
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
            #------------------------------------------------------------------------------#




            #----- Legend: the simulations. -----------------------------------------------#
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
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over sites.                                                         #
            #------------------------------------------------------------------------------#
            add = FALSE
            for (p in sequence(nsites)){
               iata = sites[p]

               #----- Skip the site if there is no data. ----------------------------------#
               ok.iata = length(obs.diel[[iata]]) > 0 && any(is.finite(obs.diel[[iata]]))
               ok.iata = ok.iata && ( ! is.na(ok.iata))
               if (ok.iata){
                  #------------------------------------------------------------------------#
                  #     We call skill twice for each site in case the site has two PCHs.   #
                  #------------------------------------------------------------------------#
                  mytaylor = taylor.plot( obs        = obs.diel[[iata]]
                                        , mod        = mod.diel[[iata]]
                                        , add        = add
                                        , pos.corr   = NA
                                        , pt.col     = simul$col
                                        , pt.bg      = simul$col
                                        , pt.pch     = sites.pch[p]
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
               #---------------------------------------------------------------------------#
            }#end for (p in sequence(nsites))
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
      }#end if (length(ref) > 2)
      #------------------------------------------------------------------------------------#
   }#end for (d in sequence(ndiel))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(ncompvar))
#------------------------------------------------------------------------------------------#

