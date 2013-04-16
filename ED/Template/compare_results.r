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
ibackground    = 2                              # Make figures compatible to background
                                                # 0 -- white
                                                # 1 -- black
                                                # 2 -- dark grey
#----- Output directory -------------------------------------------------------------------#
outroot = file.path(here,paste("structure_comp_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Site settings:                                                                       #
# eort      -- first letter ("e" or "t")                                                   #
# sites     -- site codes ("IATA")                                                         #
# sites.pch -- site symbols                                                                #
#------------------------------------------------------------------------------------------#
eort           = "e"
sites          = c("gyf","s67","s83","pdg","pnz","ban","rja","m34","cax")
sites.pch      = c(    2,    5,    9,   13,    4,    8,    1,    6,    0)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.struct     = list( name     = c("ble_iage25_pft02","ble_iage25_pft05"
                                   ,"sas_iage01_pft02","sas_iage01_pft05"
                                   ,"sas_iage25_pft02","sas_iage25_pft05"
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
outform        = c("eps","png","pdf")  # Formats for output file.  Supported formats are:
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
compvar       = list()
compvar[[ 1]] = list( vnam       = "ustar"
                    , symbol     = expression(u^symbol("\052"))
                    , desc       = "Friction velocity"
                    , unit       = untab$mos
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(purple.bg,purple.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 2]] = list( vnam       = "cflxca"
                    , symbol     = expression(F(CO[2]))
                    , desc       = "Carbon dioxide flux"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(green.bg,green.fg)
                    , leg.corner = "bottomright"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 3]] = list( vnam       = "cflxst"
                    , symbol     = expression(S(CO[2]))
                    , desc       = "Carbon dioxide storage"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(orange.bg,orange.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 4]] = list( vnam       = "nee"
                    , symbol     = expression(NEE)
                    , desc       = "Net ecosystem exchange"
                    , unit       = untab$umolcom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(green.bg,green.fg)
                    , leg.corner = "bottomright"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 5]] = list( vnam       = "nep"
                    , symbol     = expression(NEP)
                    , desc       = "Net ecosystem productivity"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(olive.bg,olive.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 6]] = list( vnam       = "reco"
                    , symbol     = expression(R[Eco])
                    , desc       = "Ecosystem respiration"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(yellow.bg,yellow.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[ 7]] = list( vnam       = "gpp"
                    , symbol     = expression(GPP)
                    , desc       = "Gross primary productivity"
                    , unit       = untab$kgcom2oyr
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(green.bg,green.fg)
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
compvar[[ 8]] = list( vnam       = "parup"
                    , symbol     = expression(PAR^symbol("\335"))
                    , desc       = "Outgoing PAR"
                    , unit       = untab$umolom2os
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(olive.bg,olive.fg)
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
compvar[[ 9]] = list( vnam       = "rshortup"
                    , symbol     = expression(SW^symbol("\335"))
                    , desc       = "Outgoing shortwave radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(indigo.bg,indigo.fg)
                    , leg.corner = "topleft"
                    , sunvar     = TRUE
                    )#end list
compvar[[10]] = list( vnam       = "rnet"
                    , symbol     = expression(R[Net])
                    , desc       = "Net radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(sky.bg,sky.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[11]] = list( vnam       = "rlongup"
                    , symbol     = expression(LW^symbol("\335"))
                    , desc       = "Outgoing longwave radiation"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(red.bg,red.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[12]] = list( vnam       = "hflxca"
                    , symbol     = expression(F(theta))
                    , desc       = "Sensible heat flux"
                    , unit       = untab$wom2
                    , col.obser  = c(grey.bg,grey.fg)
                    , col.model  = c(orange.bg,orange.fg)
                    , leg.corner = "topleft"
                    , sunvar     = FALSE
                    )#end list
compvar[[13]] = list( vnam       = "wflxca"
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
good[[ 8]] = list( vnam      = "sn.lnlike"
                 , desc      = "Scaled support based on skew normal distribution"
                 , spider    = FALSE
                 , normalise = FALSE
                 )#end list
good[[ 9]] = list( vnam      = "norm.lnlike"
                 , desc      = "Scaled support based on normal distribution"
                 , spider    = FALSE
                 , normalise = FALSE
                 )#end list
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
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
simcex.key     = rep(2.0    ,times=n.sim)
simlwd.key     = rep(2.0    ,times=n.sim)
simpch.key     = rep(21     ,times=n.sim)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Dump the information to a list.                                                      #
#------------------------------------------------------------------------------------------#
simul       = data.frame( name             = simul.key
                        , desc             = simleg.key
                        , colour           = simcol.key
                        , lty              = simlty.key
                        , cex              = simcex.key
                        , lwd              = simlwd.key
                        , pch              = simpch.key
                        , stringsAsFactors = FALSE
                        )#end data.frame
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      List the keys for all dimensions.                                                   #
#------------------------------------------------------------------------------------------#
sites.key    = sites
sites.desc   = poilist$longname[match(sites.key,poilist$iata)]
control.key  = apply(X = sapply(X=control,FUN=c),MARGIN=1,FUN=unlist)[,"vnam"]
compvar.key  = apply(X = sapply(X=compvar,FUN=c),MARGIN=1,FUN=unlist)$vnam
compvar.sym  = apply(X = sapply(X=compvar,FUN=c),MARGIN=1,FUN=unlist)$symbol
good.key     = apply(X = sapply(X=good   ,FUN=c),MARGIN=1,FUN=unlist)[,"vnam"]
season.key   = season.list
diel.key     = c("night","rise.set","day","all.hrs")
diel.desc    = c("Nighttime","Sun Rise/Set","Daytime","All hours")
#------------------------------------------------------------------------------------------#



#----- Set the various dimensions associated with variables, simulations, and sites. ------#
nsites   = length(sites.key  )
nsimul   = length(simul.key  )
ncompvar = length(compvar.key)
ncontrol = length(control.key)
ngood    = length(good.key   )
nseason  = length(season.key )
ndiel    = length(diel.key   )
#------------------------------------------------------------------------------------------#



#----- Load observations. -----------------------------------------------------------------#
obser.file = paste(srcdir,"LBA_MIP.nogapfill.RData",sep="/")
load(file=obser.file)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#----- Find the best set up for plotting all seasons in the same plot. --------------------#
lo.box = pretty.box(n=nseason-1)
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
   #     Create paths for all variables.                                                   #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ncompvar){
      this.vnam      = compvar[[v]]$vnam
      o.vnam         = list()
      o.vnam$main    = file.path(o.form$main,this.vnam)
      if (is.figure && ! file.exists(o.vnam$main)) dir.create(o.vnam$main)


      #------------------------------------------------------------------------------------#
      #     Create paths for all variables.                                                #
      #------------------------------------------------------------------------------------#
      for (d in 1:ndiel){
         this.diel      = diel.key [d]

         o.diel         = list()
         o.diel$main    = file.path(o.vnam$main,this.diel)
         o.diel$barplot = file.path(o.diel$main,"barplot")
         o.diel$skill   = file.path(o.diel$main,"skill"  )
         o.diel$taylor  = file.path(o.diel$main,"taylor" )
         if (is.figure){
            if (! file.exists(o.diel$main   )) dir.create(o.diel$main   )
            if (! file.exists(o.diel$barplot)) dir.create(o.diel$barplot)
            if (! file.exists(o.diel$skill  )) dir.create(o.diel$skill  )
            if (! file.exists(o.diel$taylor )) dir.create(o.diel$taylor )
         }#end if (is.figure)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Quality plot, we must create directories by controlling variable.            #
         #---------------------------------------------------------------------------------#
         o.qual      = list()
         o.qual$main = file.path(o.diel$main,"quality")
         if (is.figure && ! file.exists(o.qual$main)) dir.create(o.qual$main)
         for (u in 1:ncontrol){
            this.qual = control[[u]]$vnam
            o.qual[[this.qual]] = file.path(o.qual$main,this.qual)
            if (is.figure && ! file.exists(o.qual[[this.qual]])){
               dir.create(o.qual[[this.qual]])
            }#end if (is.figure && ! file.exists(o.qual[[this.qual]]))
            #------------------------------------------------------------------------------#
         }#end for (u in 1:ncontrol)
         o.diel$quality = o.qual
         #---------------------------------------------------------------------------------#



         #----- Copy the diel list to the parent list. ------------------------------------#
         o.vnam[[this.diel]] = o.diel
         #---------------------------------------------------------------------------------#
      }#end for (d in 1:ndiel)
      #------------------------------------------------------------------------------------#



      o.form[[this.vnam]]              = o.vnam
      #------------------------------------------------------------------------------------#
   }#end for for (v in 1:ncompvar)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the "spider web" plots.                                          #
   #---------------------------------------------------------------------------------------#
   o.spider                = list()
   o.spider$main           = file.path(o.form$main,"spider")
   if (is.figure && ! file.exists(o.spider$main)) dir.create(o.spider$main)
   for (d in 1:ndiel){
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


   #----- Save the full list to the main path list. ---------------------------------------#
   out[[this.form]] = o.form
   #---------------------------------------------------------------------------------------#
}#end for (o in 1:nout)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#   Loop through the sites.                                                                #
#------------------------------------------------------------------------------------------#
cat (" + Add season and diel keys to seasons and diel...","\n")
for (p in 1:nsites){
   #----- Grab the observation. -----------------------------------------------------------#
   obs        = get(paste("obs",sites[p],sep="."))
   #---------------------------------------------------------------------------------------#


   #----- Create some variables to describe season and time of the day. -------------------#
   if (! "season" %in% names(obs)) obs$season = season(obs$when,add.year=FALSE)
   if (! "diel" %in% names(obs))   obs$diel   = (! obs$nighttime) + obs$highsun
   #---------------------------------------------------------------------------------------#



   #----- Save the variables to the observations. -----------------------------------------#
   dummy = assign(paste("obs",sites[p],sep="."),obs)
   #---------------------------------------------------------------------------------------#
}#end for (p in 1:nsites)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Retrieve all data.                                                                  #
#------------------------------------------------------------------------------------------#
cat (" + Retrieve model results for all sites...","\n")
res = list()
for (p in 1:nsites){
   #----- Get the basic information. ------------------------------------------------------#
   iata          = sites[p]
   im            = match(iata,poilist$iata)
   this          = list()
   this$short    = poilist$short   [im]
   this$longname = poilist$longname[im]
   this$iata     = poilist$iata    [im]
   this$lon      = poilist$lon     [im]
   this$lat      = poilist$lat     [im]
   this$sim      = list()
   this$ans      = list()
   cat("   - Site :",this$longname,"...","\n")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Get all the statistics and actual values for every simulation.                    #
   #---------------------------------------------------------------------------------------#
   for (s in 1:nsimul){
      cat("    * Simulation: ",simul$desc[s],"...","\n")
      sim.name = paste(eort,iata,"_",simul$name[s],sep="")
      sim.path = paste(here,sim.name,sep="/")
      sim.file = paste(sim.path,"rdata_hour",paste("comp-",sim.name,".RData",sep="")
                      ,sep="/")
      load(sim.file)

      ans.name = paste("t",iata,"_",simul$name[s],sep="")
      ans.path = paste(here,sim.name,sep="/")
      ans.file = paste(sim.path,"rdata_hour",paste(sim.name,".RData",sep="")
                      ,sep="/")
      load(ans.file)
      this$sim[[simul$name[s]]] = dist.comp
      this$ans[[simul$name[s]]] = model
      rm(list=c("dist.comp","model","eddy.complete","eddy.tresume"))
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Copy the data to the results. ---------------------------------------------------#
   res[[iata]] = this
   rm(this)
   #---------------------------------------------------------------------------------------#
}#end for
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
   web = array( dim      = c   (nsimul,nsites,ncompvar,ndiel,nseason)
              , dimnames = list(simul.key,sites.key,compvar.key,diel.key,season.key)
              )#end array 
   for (v in 1:ncompvar){
      this.vnam     = compvar[[v]]$vnam
      this.measured = paste("measured",this.vnam,sep=".")

      #------------------------------------------------------------------------------------#
      #     Loop over all sites.                                                           #
      #------------------------------------------------------------------------------------#
      for (p in 1:nsites){
         iata = sites.key[p]
         sfac = matrix(data=1,nrow=ndiel,ncol=nseason,dimnames=list(diel.key,season.key))

         #---------------------------------------------------------------------------------#
         #     Find the scale factors for variables that have units.                       #
         #---------------------------------------------------------------------------------#
         if (norm.good){
            obs  = get(paste("obs",iata,sep="."))

            #----- Find out when this output variable is finite and measured. -------------#
            p.sel = is.finite(obs[[this.vnam]]) & obs[[this.measured]]
            #------------------------------------------------------------------------------#

            dd = sequence(ndiel-1)
            ee = sequence(nseason-1)

            #----- Find the components. ---------------------------------------------------#
            for (dd in sequence(ndiel)){
               d.sel = obs$diel == dd | dd == ndiel
               for (ee in sequence(nseason)){
                  e.sel = obs$season == ee | ee == nseason
                  sel   = p.sel & d.sel & e.sel
                  if (any(sel)) sfac[dd,ee] = sd(obs[[this.vnam]][sel],na.rm=TRUE)
               }#end for
            }#end for
            #------------------------------------------------------------------------------#
         }#end if (norm.good)
         sfac = ifelse(sfac == 0.,NA,1/sfac)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Grab the data for this simulation.                                          #
         #---------------------------------------------------------------------------------#
         for (s in 1:nsimul){
            this         = res[[iata]]$sim[[s]][[this.vnam]][[this.good]]
            use.season   = paste(sprintf("%2.2i",sequence(nseason)),season.key,sep="-")
            web[s,p,v,,] = abs(this[diel.key,use.season]) * sfac
         }#end for (s in 1:nsimul)
         #---------------------------------------------------------------------------------#
      }#end for (p in 1:nsites)
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:ncompvar)
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Plot the spider webs by diel and variable.                                        #
   #---------------------------------------------------------------------------------------#
   for (d in 1:ndiel){
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Webs by site (all variables).                                                  #
      #------------------------------------------------------------------------------------#
      for (p in 1:nsites){
         iata     = sites.key[p]

         letitre = paste(desc.good," - ",sites.desc[p],"\n",diel.desc[d],sep="")

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

            #------------------------------------------------------------------------------#
            #     Webs by variable (all sites).                                            #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Make the file name. -------------------------------------------------#
               out.web = out[[outform[o]]]$spider[[diel.key[d]]]$sites
               fichier   = paste(out.web,"/spider-allyear-",iata,"-",this.good
                                        ,"-",diel.key[d],".",outform[o],sep="")
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
                      , cex     = cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the spider web.                                                  #
               #---------------------------------------------------------------------------#
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
      #     Webs by variable (all sites).                                                  #
      #------------------------------------------------------------------------------------#
      for (v in 1:ncompvar){
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

            #------------------------------------------------------------------------------#
            #     Webs by variable (all sites).                                            #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Make the file name. -------------------------------------------------#
               out.web = out[[outform[o]]]$spider[[diel.key[d]]]$variables
               fichier   = paste(out.web,"/spider-allyear-",this.vnam,"-",this.good
                                        ,"-",diel.key[d],".",outform[o],sep="")
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
                      , pt.cex  = simul$cex
                      , cex     = cex.ptsz
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the spider web.                                                  #
               #---------------------------------------------------------------------------#
               radial.flex( lengths          = web[,p.sel,v,d,nseason]
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

            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end if (any(is.finite(web[,,v,d,nseason])))
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:ncompvar)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



   }#end for (d in 1:ndiel)
   #---------------------------------------------------------------------------------------#



}#end for (g in good.loop)
#------------------------------------------------------------------------------------------#







#------------------------------------------------------------------------------------------#
#         Plot the various statistics as functions of the site "completion".               #
#------------------------------------------------------------------------------------------#
cat (" + Plot statistics as functions of model and fraction of input data...","\n")
performance = array( data     = NA
                   , dim      = c(ncompvar,ngood,nsimul)
                   , dimnames = list(compvar.key,good.key,simul.key)
                   )#end array
for (v in 1:ncompvar){
   #----- Copy the variable information. --------------------------------------------------#
   this.vnam     = compvar[[v]]$vnam
   this.desc     = compvar[[v]]$desc
   this.unit     = compvar[[v]]$unit
   this.sun      = compvar[[v]]$sunvar
   this.measured = paste("measured",this.vnam,sep=".")
   cat("   - ",this.desc,"...","\n")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Loop over all sites, seasons, and diel and get the average score for all input    #
   # variables for when the observations are valid.                                        #
   #---------------------------------------------------------------------------------------#
   input.score = array( data     = 0.
                      , dim      = c(ndiel,nseason,ncontrol,nsites)
                      , dimnames = list(diel.key,season.key,control.key,sites.key)
                      )#end score
   output.nobs = array( data     = 0.
                      , dim      = c(ndiel,nseason,nsites)
                      , dimnames = list(diel.key,season.key,sites.key)
                      )#end score
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Loop through the sites.                                                             #
   #---------------------------------------------------------------------------------------#
   for (p in 1:nsites){
      #----- Grab the observation. --------------------------------------------------------#
      obs        = get(paste("obs",sites[p],sep="."))
      #------------------------------------------------------------------------------------#


      #----- Find out when this output variable is finite and measured. -------------------#
      p.sel = is.finite(obs[[this.vnam]]) & obs[[this.measured]]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Loop over all seasons.                                                        #
      #------------------------------------------------------------------------------------#
      for (e in 1:nseason){
         #----- Select this season (or everything for all seasons). -----------------------#
         e.sel = obs$season == e | e == nseason
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all parts of the day.                                            #
         #---------------------------------------------------------------------------------#
         for (d in 1:ndiel){
            #----- Select this diel (or everything for all day). --------------------------#
            d.sel = obs$diel == (d-1) | d == ndiel
            d.sel = d.sel & ((! this.sun) | obs$highsun)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Combine the selections.                                                 #
            #------------------------------------------------------------------------------#
            sel   = e.sel & d.sel & p.sel
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Loop over all control variables (except the "global" one).              #
            #------------------------------------------------------------------------------#
            for (u in 1:(ncontrol-1)){
               #---- Get the score for this variable. -------------------------------------#
               control.vnam = paste("score",control[[u]]$vnam,sep=".")
               obs.score    = obs[[control.vnam]]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Combine all the selections and find the mean score.                  #
               #---------------------------------------------------------------------------#
               input.score[d,e,u,p] = mean(obs.score[sel],na.rm=TRUE)
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Find the general score and number of observations.                      #
            #------------------------------------------------------------------------------#
            input.score[d,e,ncontrol,p] = mean(input.score[d,e,1:(ncontrol-1),p],na.rm=TRUE)
            output.nobs[d,e,p]          = sum(sel,na.rm=TRUE)
            #------------------------------------------------------------------------------#
         }#end for (d in 1:ndiel)
         #---------------------------------------------------------------------------------#
      }#end for (e in 1:nseason)
      #------------------------------------------------------------------------------------#
   }#end for (p in 1:nsites)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make sure that non-finite scores are NA.                                         #
   #---------------------------------------------------------------------------------------#
   input.score[! is.finite(input.score)] = NA
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Loop over all parts of the day.                                                  #
   #---------------------------------------------------------------------------------------#
   cat("     * Skill and Taylor plots...","\n")
   for (d in 1:ndiel){
      cat("       > ",diel.desc[d],"...","\n")


      #------------------------------------------------------------------------------------#
      #      Loop over all sites, normalise the data and create the vector for the model.  #
      #------------------------------------------------------------------------------------#
      obs.diel    = list()
      mod.diel    = list()
      bias.range  = NULL
      sigma.range = NULL
      for (p in 1:nsites){
         iata  = sites[p]
         obs   = get(paste("obs",iata,sep="."))
         nwhen = length(obs$when)


         #----- Select this diel (or everything for all day). -----------------------------#
         d.sel = (obs$diel == (d-1) | d == ndiel) & ((! this.sun) | obs$highsun)
         sel   = d.sel & is.finite(obs[[this.vnam]]) & obs[[this.measured]]
         sel   = sel   & is.finite(sel)
         n.sel = sum(sel)
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Find the standard deviation of this observation.  Skip the site if every-   #
         # thing is zero.                                                                  #
         #---------------------------------------------------------------------------------#
         this.obs     = obs[[this.vnam]][sel]
         sdev.obs.now = sd(this.obs,na.rm=TRUE)
         sel          = sel & is.finite(sdev.obs.now) & sdev.obs.now > 0
         #---------------------------------------------------------------------------------#



         #----- Copy the observed data. ---------------------------------------------------#
         obs.diel[[iata]] = this.obs
         #---------------------------------------------------------------------------------#



         #----- Copy the modelled data, and update ranges. --------------------------------#
         mod.diel[[iata]] = matrix(ncol=nsimul,nrow=n.sel,dimnames=list(NULL,simul.key))
         if (any(sel)){
            for (s in 1:nsimul){
               this.mod             = res[[iata]]$ans[[simul.key[s]]][[this.vnam]][sel]
               this.res             = this.mod - this.obs
               mod.diel[[iata]][,s] = this.mod

               #----- Find the normalised bias and model standard deviation. --------------#
               bias.now    = mean(this.res, na.rm=TRUE) / sdev.obs.now
               sigma.now   = sd  (this.res, na.rm=TRUE) / sdev.obs.now
               bias.range  = c(bias.range ,bias.now   )
               sigma.range = c(sigma.range,sigma.now  )
               if (! is.finite(bias.now) || ! is.finite(sigma.now)){
                  stop("Something weird with this run...")
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (s in 1:nsimul)
            #------------------------------------------------------------------------------#
         }#end if (any(sel))
         #---------------------------------------------------------------------------------#
      }#end for (p in 1:nsites)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Plot Taylor and skill plots only if there is anything to plot.                 #
      #------------------------------------------------------------------------------------#
      if (length(unlist(obs.diel)) > 0){
         #---- Fix ranges. ----------------------------------------------------------------#
         xy.range    = 1.04 * max(c(abs(bias.range),abs(sigma.range)))
         bias.range  = 1.04 * xy.range  * c(-1,1)
         sigma.range = 1.04 * xy.range  * c( 1,0)
         r2.range    = range(1-xy.range^2,1)
         #---------------------------------------------------------------------------------#



         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Skill plot.                                                                 #
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Plot title.                                                                 #
         #---------------------------------------------------------------------------------#
         letitre = paste(" Skill diagram - ",this.desc,"\n",diel.desc[d],sep="")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Make the file name. ----------------------------------------------------#
            out.skill = out[[outform[o]]][[this.vnam]][[diel.key[d]]]$skill
            fichier   = paste(out.skill,"/skill-allyear-",this.vnam,"-",diel.key[d]
                             ,".",outform[o],sep="")
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
            layout(mat = rbind(c(3,3,3),c(1,2,2)),height = c(5.0,1.0))
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
                   , pt.cex  = mean(unique(simul$cex))
                   , cex     = 1.1 * cex.ptsz
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
                   , pch     = simul$pch
                   , ncol    = 2
                   , title   = expression(bold("Structure"))
                   , cex     = 0.9 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over sites.                                                         #
            #------------------------------------------------------------------------------#
            myskill = NULL
            for (p in 1:nsites){
               iata = sites[p]


               #----- Skip the site if there is no data. ----------------------------------#
               if (length(obs.diel[[iata]]) > 0){
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
                                                            , cex = simul$cex
                                                            , lty = "solid"
                                                            , lwd = simul$lwd
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
         for (o in 1:nout){
            #----- Make the file name. ----------------------------------------------------#
            out.taylor = out[[outform[o]]][[this.vnam]][[diel.key[d]]]$taylor
            fichier    = paste(out.taylor,"/taylor-allyear-",this.vnam,"-",diel.key[d]
                              ,".",outform[o],sep="")
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
            layout(mat = rbind(c(3,3,3),c(1,2,2)),height = c(5.0,1.0))
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
                   , pt.cex  = 2.0
                   , cex     = cex.ptsz
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
                   , pch     = simul$pch
                   , ncol    = 2
                   , title   = expression(bold("Structure"))
                   , cex     = 14/ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over sites.                                                         #
            #------------------------------------------------------------------------------#
            add = FALSE
            for (p in 1:nsites){
               iata = sites[p]

               #----- Skip the site if there is no data. ----------------------------------#
               if (length(obs.diel[[iata]]) > 0.){
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
                                        , pt.cex     = simul$cex
                                        , pt.lwd     = simul$lwd
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
      }#end if (length(ref) > 2)
      #------------------------------------------------------------------------------------#
   }#end for (d in 1:ndiel)
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      Loop over all sites and simulations, and grab the data.                          #
   #---------------------------------------------------------------------------------------#
   for (g in 1:ngood){
      nobs.good = output.nobs


      this.good = good[[g]]$vnam
      desc.good = good[[g]]$desc
      cat ("     * ",desc.good,"\n")
      stat    = array( data     = NA
                     , dim      = c(ndiel,nseason,nsimul,nsites)
                     , dimnames = list(diel.key,season.key,simul.key,sites.key)
                     )#end array
      nvars   = array( data     = NA
                     , dim      = c(ndiel,nseason,nsites)
                     , dimnames = list(diel.key,season.key,sites.key)
                     )#end array
      colstat = array( data     = NA_character_
                     , dim      = c(nsimul,nsites)
                     , dimnames = list(simul.key,sites.key)
                     )#end array
      #------------------------------------------------------------------------------------#
      #      Loop over all sites and simulations.                                          #
      #------------------------------------------------------------------------------------#
      for (p in 1:nsites){
         iata         = sites[p]
         obs          = get(paste("obs",sites[p],sep="."))
         for (s in 1:nsimul){
            colstat[s,p]   = simul$colour[s]
            this           = res[[iata]]$sim[[s]][[this.vnam]][[this.good]]

            use.season  = paste(sprintf("%2.2i",sequence(nseason)),season.key,sep="-")
            stat[,,s,p] = this[diel.key,use.season]
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      bye            = apply(X=(! is.finite(stat)),MARGIN=c(1,2,4),FUN=sum,na.rm=TRUE)
      bye            = bye != 0
      nobs.good[bye] = -1
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Standardise the log-likelihood so the data are more comparable.                #
      #------------------------------------------------------------------------------------#
      if (this.good %in% c("lsq.lnlike","sn.lnlike","norm.lnlike")){
         stat.min  = apply(X = stat, MARGIN=c(1,2,4),FUN=min,na.rm=TRUE)
         stat.max  = apply(X = stat, MARGIN=c(1,2,4),FUN=max,na.rm=TRUE)
         stat.orig = stat
         for (s in 1:nsimul){
           stat[,,s,] = 100. * ( (stat.orig[,,s,] - stat.max)/ (stat.max - stat.min) )
         }#end for
      }else if(this.good %in% "r.squared"){
         orig.stat  = stat
         sel        = is.finite(stat) & stat < -1
         stat[sel]  = -1
      }#end if
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #       Plot the box plots.                                                          #
      #------------------------------------------------------------------------------------#
      for (d in 1:ndiel){
         #---------------------------------------------------------------------------------#
         #      Find the limits for the bar plot.                                          #
         #---------------------------------------------------------------------------------#
         xlimit  = c(0,nsites*(nsimul+1))+0.5
         xat     = seq(from=0,to=(nsites-1)*(nsimul+1),by=nsimul+1)+1+0.5*nsimul
         xlines  = seq(from=0,to=nsites*(nsimul+1),by=nsimul+1)+0.5
         if (this.good %in% c("r.squared")){
            y.nobs =  1.10
            y.r2   = -1.10
            ylimit = c(-1.2,1.2)
         }else{
            ylimit  = range(stat[d,,,],na.rm=TRUE)
            if (  any(! is.finite(ylimit)) || (ylimit[1] == ylimit[2] && ylimit[1] == 0)){
               y.nobs    = 0.90
               ylimit    = c(-1,1)
            }else if (ylimit[1] == ylimit[2] ){
               y.nobs    = ylimit[2] * ( 1. + sign(ylimit[2]) * 0.50 * fracexp )
               ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1])        * fracexp )
               ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2])        * fracexp )
            }else{
               y.nobs    = ylimit[2] + 0.15 * fracexp * diff(ylimit)
               ylimit[2] = ylimit[2] + 0.30 * fracexp * diff(ylimit)
            }#end if
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Loop over all output formats.                                               #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Make the file name. ----------------------------------------------------#
            out.barplot = out[[outform[o]]][[this.vnam]][[diel.key[d]]]$barplot
            fichier     = paste(out.barplot,"/bplot-byseason-",this.vnam,"-",this.good,"-"
                                           ,diel.key[d],".",outform[o],sep="")
            if (outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth
                  ,height=size$height*depth,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,4,0))
            layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                  ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simleg.key
                   , fill    = simcol.key
                   , border  = foreground
                   , ncol    = min(3,pretty.box(nsimul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all seasons, and plot the bar plots.                           #
            #------------------------------------------------------------------------------#
            for (e in 1:(nseason-1)){
               #----- Find out where is this box going, and set up axes and margins. ------#
               left    = (e %% lo.box$ncol) == 1
               right   = (e %% lo.box$ncol) == 0
               top     = e <= lo.box$ncol
               bottom  = e > (lo.box$nrow - 1) * lo.box$ncol
               mar.now = c(4 + 0 * bottom,1 + 1 * left,3 + 0 * top,1 + 1 * right) + 0.1
               #---------------------------------------------------------------------------#


               #----- Set up the title for each plot. -------------------------------------#
               lesub = paste(season.full[e],sep="")
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Order the sites by amount of information.                             #
               #---------------------------------------------------------------------------#
               op   = order(nobs.good[d,e,],na.last=FALSE)
               #---------------------------------------------------------------------------#



               #----- Plot window and grid. -----------------------------------------------#
               par(mar=mar.now,xpd=FALSE,las=2)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
               axis(side=1,at=xat,labels=toupper(sites.key[op]))
               if (left  ) axis(side=2)
               box()
               title(main=lesub)
               if (plotgrid) abline(h=axTicks(2),v=xlines,col=grid.colour,lty="solid")
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Add the bar plot.                                                     #
               #---------------------------------------------------------------------------#
               xbp = barplot(height=stat[d,e,,op],col=colstat,beside=TRUE,border=grey.fg
                            ,add=TRUE,axes=FALSE,axisnames=FALSE,xpd=FALSE)
               text   (x=xat,y=y.nobs,labels=nobs.good[d,e,op],cex=0.7)
               if (this.good %in% c("r.squared")){
                  ybp = y.r2 - 100 * ( orig.stat[d,e,,op] >= -1.0)
                  text(x=xbp,y=ybp,labels=rep("*",times=length(ybp))
                      ,font=2,col="sandybrown")
               }#end if
               #---------------------------------------------------------------------------#
            }#end for



            #------------------------------------------------------------------------------#
            #     Make the title and axis labels.                                          #
            #------------------------------------------------------------------------------#
            letitre = paste(desc.good," - ",this.desc,"\n",diel.desc[d],sep="")
            if (this.good %in% c("bias","rmse")){
               ley  = desc.unit(desc=desc.good,unit=this.unit)
            }else{
               ley  = desc.unit(desc=desc.good,unit=untab$empty)
            }#end if
            lex     = "Sites"
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=ley    ,side=2,outer=TRUE,padj=-0.75)
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

         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot all seasons together.                                                  #
         #---------------------------------------------------------------------------------#
         #----- Make the title and axis labels. -------------------------------------------#
         letitre = paste(desc.good," - ",this.desc,"\n",diel.desc[d]," - All seasons"
                        ,sep="")
         if (this.good %in% c("bias","rmse")){
            ley  = desc.unit(desc=desc.good,unit=this.unit)
         }else{
            ley  = desc.unit(desc=desc.good,unit=untab$empty)
         }#end if
         lex     = "Sites"
         #-----  Order the sites by amount of information. --------------------------------#
         nobs = nobs.good[d,nseason,]
         op   = order(nobs,na.last=FALSE)
         #----- Loop over all formats. ----------------------------------------------------#
         for (o in 1:nout){
            #----- Make the file name. ----------------------------------------------------#
            out.barplot = out[[outform[o]]][[this.vnam]][[diel.key[d]]]$barplot
            fichier     = paste(out.barplot,"/bplot-allyear-",this.vnam,"-",this.good,"-"
                                           ,diel.key[d],".",outform[o],sep="")
            if (outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth
                  ,height=size$height*depth,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=size$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE
                  ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            layout(mat = rbind(2,1), height = c(5,1))
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,4.1,0.1,4.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "top"
                   , inset   = 0.0
                   , legend  = simleg.key
                   , fill    = simcol.key
                   , border  = foreground
                   , ncol    = min(3,pretty.box(nsimul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 0.85 * cex.ptsz
                   , xpd     = TRUE
                   )#end legend
            #------------------------------------------------------------------------------#


            #----- Set up the title for each plot. ----------------------------------------#
            lesub = paste("All seasons",sep="")
            #------------------------------------------------------------------------------#



            #----- Plot window and grid. --------------------------------------------------#
            par(mar=c(5,4,3,2)+0.1,xpd=FALSE,las=2)
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
            title(main=letitre,xlab=lex,ylab=ley)
            axis(side=1,at=xat,labels=toupper(sites.key[op]))
            axis(side=2)
            box()
            if (plotgrid) abline(h=axTicks(2),v=xlines,col=grid.colour,lty="solid")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Add the bar plot.                                                        #
            #------------------------------------------------------------------------------#
            barplot(height=stat[d,nseason,,op],col=colstat,beside=TRUE,border=grey.fg
                   ,add=TRUE,axes=FALSE,axisnames=FALSE,xpd=FALSE)
            text(x=xat,y=y.nobs,labels=nobs.good[d,nseason,op],cex=0.9,col=grey.fg)
            if (this.good %in% c("r.squared")){
               ybp = y.r2 - 100 * ( orig.stat[d,nseason,,op] >= -1.0)
               text(x=xbp,y=ybp,labels=rep("*",times=length(ybp)),font=2,col="sandybrown")
            }#end if
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

         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Plot the statistics as a function of the quality of the data.              #
         #---------------------------------------------------------------------------------#
         for (u in 1:ncontrol){
            this.qual = control[[u]]
            qual.vnam = this.qual$vnam
            qual.desc = this.qual$desc
            qual.unit = this.qual$unit



            #------------------------------------------------------------------------------#
            #      Find the limits for the bar plot.                                       #
            #------------------------------------------------------------------------------#
            xlimit  = range(input.score[d,,u,],na.rm=TRUE)
            ylimit  = range(stat       [d,,, ],na.rm=TRUE)
            if (  any(! is.finite(xlimit)) || (xlimit[1] == xlimit[2] && xlimit[1] == 0)){
               xlimit    = c(0,10)
            }else if (xlimit[1] == xlimit[2] ){
               xlimit[1] = xlimit[1] * ( 1. - sign(xlimit[1]) * fracexp)
               xlimit[2] = xlimit[2] * ( 1. + sign(xlimit[2]) * fracexp)
            }#end if
            if (this.good %in% c("r.squared")){
               ylimit = c(-1,1)
            }else{
               ylimit  = range(stat[d,,,],na.rm=TRUE)
               if (  any(! is.finite(ylimit)) 
                  || (ylimit[1] == ylimit[2] && ylimit[1] == 0) ){
                  ylimit    = c(-1,1)
               }else if (ylimit[1] == ylimit[2] ){
                  ylimit[1] = ylimit[1] * ( 1. - sign(ylimit[1]) * fracexp)
                  ylimit[2] = ylimit[2] * ( 1. + sign(ylimit[2]) * fracexp)
               }#end if
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all output formats.                                            #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Make the file name. -------------------------------------------------#
               out.diel    = out[[outform[o]]][[this.vnam]][[diel.key[d]]]
               out.qualnow = out.diel$quality[[qual.vnam]]
               fichier     = paste(out.qualnow,"/qual-byseason-",this.vnam,"-",qual.vnam
                                  ,"-",this.good,"-",diel.key[d],".",outform[o],sep="")
               if (outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth
                     ,height=size$height*depth,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into several smaller windows.  Add a bottom row to   #
               # fit the legend.                                                           #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,4,0))
               layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                     ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#



               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = simleg.key
                      , fill    = simcol.key
                      , border  = foreground
                      , ncol    = min(3,pretty.box(nsimul)$ncol)
                      , title   = expression(bold("Simulation"))
                      , cex     = cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over all seasons, and plot the bar plots.                        #
               #---------------------------------------------------------------------------#
               for (e in 1:(nseason-1)){
                  #----- Find out where is this box going, and set up axes and margins. ---#
                  left    = (e %% lo.box$ncol) == 1
                  right   = (e %% lo.box$ncol) == 0
                  top     = e <= lo.box$ncol
                  bottom  = e > (lo.box$nrow - 1) * lo.box$ncol
                  mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
                  #------------------------------------------------------------------------#


                  #----- Set up the title for each plot. ----------------------------------#
                  lesub = paste(season.full[e],sep="")
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=mar.now,xpd=FALSE,las=1)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  if (bottom) axis(side=1)
                  if (left  ) axis(side=2)
                  box()
                  title(main=lesub)
                  if (plotgrid) grid(col=grid.colour,lty="solid")
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Add the bar plot.                                                  #
                  #------------------------------------------------------------------------#
                  for (s in 1:nsimul){
                     points(x=input.score[d,e,u,],y=stat[d,e,s,],pch=15
                           ,col=simul$colour[s],cex=1.0)
                  }#end for
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Make the title and axis labels.                                       #
               #---------------------------------------------------------------------------#
               letitre = paste(desc.good," - ",this.desc,"\n",diel.desc[d],sep="")
               if (this.good %in% c("bias","rmse")){
                  ley  = desc.unit(desc=desc.good,unit=this.unit)
               }else{
                  ley  = desc.unit(desc=desc.good,unit=untab$empty)
               }#end if
               lex     = paste("Quality index - ",qual.desc,sep="")
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=lex    ,side=1,outer=TRUE,padj=-6.00)
               mtext(text=ley    ,side=2,outer=TRUE,padj=-0.75)
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






            #------------------------------------------------------------------------------#
            #      Plot the combined data for all seasons.                                 #
            #------------------------------------------------------------------------------#
            #----- Make the title and axis labels. ----------------------------------------#
            letitre = paste(desc.good," - ",this.desc,"\n",diel.desc[d]," - All seasons"
                           ,sep="")
            if (this.good %in% c("bias","rmse")){
               ley  = desc.unit(desc=desc.good,unit=this.unit)
            }else{
               ley  = desc.unit(desc=desc.good,unit=untab$empty)
            }#end if
            lex     = paste("Quality index - ",qual.desc,sep="")
            #----- Loop over all output formats. ------------------------------------------#
            for (o in 1:nout){
               #----- Make the file name. -------------------------------------------------#
               out.diel    = out[[outform[o]]][[this.vnam]][[diel.key[d]]]
               out.qualnow = out.diel$quality[[qual.vnam]]
               fichier     = paste(out.qualnow,"/qual-allyear-",this.vnam,"-",qual.vnam
                                  ,"-",this.good,"-",diel.key[d],".",outform[o],sep="")
               if (outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth
                     ,height=size$height*depth,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into several smaller windows.  Add a bottom row to   #
               # fit the legend.                                                           #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               layout(mat=rbind(2,1),height=c(5,1))
               #---------------------------------------------------------------------------#



               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,4.1,0.1,4.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "top"
                      , inset   = 0.0
                      , legend  = simleg.key
                      , fill    = simcol.key
                      , border  = foreground
                      , ncol    = min(3,pretty.box(nsimul)$ncol)
                      , title   = expression(bold("Simulation"))
                      , cex     = cex.ptsz
                      , xpd     = TRUE
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the annual data, with not distinction between day and night.     #
               #---------------------------------------------------------------------------#
               #----- Plot window and grid. -----------------------------------------------#
               par(mar=c(5,4,3,2)+0.1,xpd=FALSE,las=1)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
               title(main=letitre,xlab=lex,ylab=ley)
               axis(side=1)
               axis(side=2)
               box()
               if (plotgrid) grid(col=grid.colour,lty="solid")
               #---------------------------------------------------------------------------#



               #---- Add the bar plot. ----------------------------------------------------#
               for (s in 1:nsimul){
                  points(x=input.score[d,nseason,u,],y=stat[d,nseason,s,],pch=15
                        ,col=simul$colour[s],cex=1.0)
               }#end for
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





         }#end for (u in 1:ncontrol)
         #---------------------------------------------------------------------------------#
      }#end for (d in 1:ndiel)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find the general score for all variables.                                     #
      #------------------------------------------------------------------------------------#
      for (s in 1:nsimul){
         performance[v,g,s] = weighted.mean( x = stat       [,,s       ,]
                                           , w = input.score[,,ncontrol,] 
                                               * nobs.good  [,,         ]
                                           , na.rm = TRUE
                                           )#end weighted.mean
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for (g in 1:ngood)
   #---------------------------------------------------------------------------------------#
}#end for (v in 1:ncompvar)
#------------------------------------------------------------------------------------------#

