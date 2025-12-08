#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
here        = getwd()                           #   Current directory
srcdir      = "/prj/prjidfca/marcosl/Util/Rsc"  #   Script directory
ibackground = 0                                 # Make figures compatible to background
                                                # 0 -- white
                                                # 1 -- black
                                                # 2 -- dark grey
#----- Output directory -------------------------------------------------------------------#
outroot = file.path(here,paste("tseries_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform           = c("pdf")  # Formats for output file.  Supported formats are:
                              #   - "X11"    - for printing on screen (Linux)
                              #   - "quartz" - printing on screen (Mac)
                              #   - "eps"    - for postscript printing
                              #   - "png"    - for PNG printing
                              #   - "pdf"    - for PDF printing
depth             = 96        # PNG resolution, in pixels per inch
wpaper            = "letter"  # Wide paper size, to define the plot shape
spaper            = "square"  # Square paper size, to define the plot shape
plotgrid          = FALSE     # Plot grid?
multiple.lty      = FALSE     # Use different line types
ptsz              = 18        # Font size.
f.leg             = 1/6       # Expansion factor to fit legend.
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Define some dimensions.                                                              #
#------------------------------------------------------------------------------------------#
pft.use = c(2, 3,4,18)                  # PFTs to include (add PFT=18, the total)
pft.lty = c("dashed","dotdash","solid") # PFT line types
#------------------------------------------------------------------------------------------#





#------ Miscellaneous settings. -----------------------------------------------------------#
yeara          = 2002         # First year we will include
yearz          = 2007         # Last year we will include
out.monthly    = TRUE         # Output monthly averages (FALSE means yearly averages)
ystep          = 1            # Time step for PFT/DBH plots.
idbh.type      = 3            # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                              # 3 -- 0-10; 10-35; 35-70; > 70 (cm)
slz.use        = -0.10        # Maximum depth to consider for soil variables
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
#                                                                                          #
# simul     -- a list of simulations that could be processed.  Each list element is also   #
#              a list with the following information:                                      #
#              name      -- simulation name                                                #
#              desc      -- simulation description (for titles and legends)                #
#              lty       -- line type (ignored if multiple.lty is FALSE)                   #
#              col       -- line colour                                                    #
#------------------------------------------------------------------------------------------#
eort  = "t"
#----- Sites. -----------------------------------------------------------------------------#
s     = 0
sites = list()
s          = s + 1
sites[[s]] = list( iata = "s67", desc = "Tapajos National Forest")
#s          = s + 1
#sites[[s]] = list( iata = "tnf", desc = "Tapajos National Forest")
#----- Simulations. -----------------------------------------------------------------------#
u          = 0
simul      = list()
u          = u + 1
simul[[u]] = list( name  = "ivdyn00_ihrzrad00"
                 , desc  = "Dynamics OFF; Horizontal OFF"
                 , lty   = "solid"
                 , col   = "#3B24B3"
                 )#end list
u          = u + 1
simul[[u]] = list( name = "ivdyn00_ihrzrad01"
                 , desc = "Dynamics OFF; Horizontal ON"
                 , lty  = "longdash"
                 , col  = "#2996CC"
                 )#end list
u          = u + 1
simul[[u]] = list( name = "ivdyn01_ihrzrad00"
                 , desc = "Dynamics ON; Horizontal OFF"
                 , lty  = "dotdash"
                 , col  = "#990F0F"
                 )#end list
u          = u + 1
simul[[u]] = list( name = "ivdyn01_ihrzrad01"
                 , desc = "Dynamics ON; Horizontal ON"
                 , lty  = "dotted"
                 , col  = "#E65C17"
                 )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Variable settings:                                                                   #
#------------------------------------------------------------------------------------------#
n     = 0
variables      = list()
n              = n + 1
variables[[n]] = list( vnam   = "fire.dist"
                     , desc   = "Fire disturbance rate"
                     , unit   = "oneoyr"
                     , dset   = "mdist"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = TRUE
                     , ilu    = 4
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "agb"
                     , desc   = "Above-ground biomass"
                     , unit   = "kgcom2"
                     , dset   = "szpft"
                     , soil   = FALSE
                     , census = TRUE
                     , zero   = TRUE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "ba"
                     , desc   = "Basal area"
                     , unit   = "cm2om2"
                     , dset   = "szpft"
                     , soil   = FALSE
                     , census = TRUE
                     , zero   = TRUE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "lai"
                     , desc   = "Leaf area index"
                     , unit   = "m2lom2"
                     , dset   = "szpft"
                     , soil   = FALSE
                     , census = TRUE
                     , zero   = TRUE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "gpp"
                     , desc   = "Gross primary productivity"
                     , unit   = "kgcom2oyr"
                     , dset   = "szpft"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "plant.resp"
                     , desc   = "Plant respiration"
                     , unit   = "kgcom2oyr"
                     , dset   = "szpft"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "agb.change"
                     , desc   = "AGB change rate"
                     , unit   = "pcagboyr"
                     , dset   = "szpft"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "hflxca"
                     , desc   = "Sensible heat flux"
                     , unit   = "Wom2"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "wflxca"
                     , desc   = "Water flux"
                     , unit   = "kgwom2oday"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "fast.soil.c"
                     , desc   = "Fast Soil Carbon"
                     , unit   = "kgcom2"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = TRUE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "struct.soil.c"
                     , desc   = "Structural Soil Carbon"
                     , unit   = "kgcom2"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = TRUE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "slow.soil.c"
                     , desc   = "Slow Soil Carbon"
                     , unit   = "kgcom2"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = TRUE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "soil.water"
                     , desc   = "Soil moisture (top 10cm)"
                     , unit   = "m3wom3"
                     , dset   = "emean"
                     , soil   = TRUE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "soil.temp"
                     , desc   = "Soil moisture (top 10cm)"
                     , unit   = "degC"
                     , dset   = "emean"
                     , soil   = TRUE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "can.temp"
                     , desc   = "Canopy air space temperature"
                     , unit   = "degC"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "rshort.gnd"
                     , desc   = "Ground absorbed SW radiation"
                     , unit   = "wom2"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "par.gnd"
                     , desc   = "Ground absorbed PAR"
                     , unit   = "umolom2os"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "par.leaf"
                     , desc   = "Absolute leaf absorption - PAR"
                     , unit   = "umolom2os"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "leaf.par"
                     , desc   = "Norm. leaf absorption - PAR"
                     , unit   = "umolom2os"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "leaf.gsw"
                     , desc   = "Stomatal conductance"
                     , unit   = "kgwom2loday"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "leaf.gpp"
                     , desc   = "Leaf-level GPP"
                     , unit   = "kgcom2loyr"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "assim.light"
                     , desc   = "Light-limited Assim. Rate"
                     , unit   = "umolom2los"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
n              = n + 1
variables[[n]] = list( vnam   = "assim.rubp"
                     , desc   = "RuBP-limited Assim. Rate"
                     , unit   = "umolom2los"
                     , dset   = "emean"
                     , soil   = FALSE
                     , census = FALSE
                     , zero   = FALSE
                     , ilu    = NA
                     )#end list
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


#----- Ensure output path exists. ---------------------------------------------------------#
if (! file.exists(outroot )) dir.create(outroot )
#------------------------------------------------------------------------------------------#


#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (1. - f.leg)
wsize   = plotsize(proje=FALSE,paper=wpaper,extendfc="lat",extfactor=f.ext)
ssize   = plotsize(proje=FALSE,extendfc="lat",extfactor=f.ext,paper=spaper)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Convert sites to data frame, and keep only those to be run.                         #
#------------------------------------------------------------------------------------------#
sites     = list.2.data.frame(sites    )
simul     = list.2.data.frame(simul    )
variables = list.2.data.frame(variables)
if (! multiple.lty) simul$lty = "solid"
#------------------------------------------------------------------------------------------#



#----- Alias to some useful pft definitions. ----------------------------------------------#
pft.key       = pft$key    [pft.use]            # PFT keys  (for dimnames)
pft.desc      = pft$name   [pft.use]            # PFT names (for titles)
pft.colour    = pft$colour [pft.use]            # PFT colours
npft.use  = length(pft.use)
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Define some dimensions.                                                              #
#------------------------------------------------------------------------------------------#
dbh.use    = sequence(ndbh+1)
dbh.key    = c(dbhkeys,"all")
dbh.desc   = c(dbhnames,"All")
dbh.colour = dbhcols
dbh.lty    = c("dotted","dotdash","solid") # PFT line types
ndbh.use   = length(dbh.use)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Convert sites to data frame, and keep only those to be run.                         #
#------------------------------------------------------------------------------------------#
nsites     = nrow(sites    )
nsimul     = nrow(simul    )
nvariables = nrow(variables)
nout       = length(outform)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Find monthly time series.                                                           #
#------------------------------------------------------------------------------------------#
nyears    = yearz-yeara+1
when      = chron( paste( rep(x=sequence(12),times=nyears)
                        , 1
                        , rep(x=seq(from=yeara,to=yearz,by=1),each=12)
                        , sep = "/"
                        )#end paste
                 )#end chron
wlabel    = as.character(when)
whenmonth = nummonths(when)
whenyear  = numyears(when)
nwhen     = length(when)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Initialise arrays                                                                    #
#------------------------------------------------------------------------------------------#
output   = array( dim      = c(nwhen,nvariables,nsimul,nsites)
                , dimnames = list(wlabel,variables$vnam,simul$name,sites$iata)
                )#end array
outsss   = array( dim      = c(nwhen,ndbh.use,npft.use,nvariables,nsimul,nsites)
                , dimnames = list(wlabel,dbh.key,pft.key,variables$vnam
                                 ,simul$name,sites$iata)
                )#end array
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Loop over sites, simulations, and variables.                                        #
#------------------------------------------------------------------------------------------#
cat(" + Load data sets...","\n")
for (s in sequence(nsites)){
   #----- Alias for 'sites' variables. ----------------------------------------------------#
   s.iata = sites$iata[s]
   s.desc = sites$desc[s]
   #---------------------------------------------------------------------------------------#

   #----- Loop over simulations. ----------------------------------------------------------#
   for (u in sequence(nsimul)){
      #----- Alias for 'simul' variables. -------------------------------------------------#
      u.name = simul$name[u]
      u.desc = simul$desc[u]
      #------------------------------------------------------------------------------------#


      #----- Load files. ------------------------------------------------------------------#
      prefix = paste(eort,s.iata,"_",u.name,sep="")
      Rfile  = file.path(here,prefix,"rdata_month",paste0(prefix,".RData"))
      cat ("   - ",s.desc," (",basename(Rfile),")","...","\n")
      load(file=Rfile)
      #------------------------------------------------------------------------------------#


      #----- Drop datum. ------------------------------------------------------------------#
      emean = datum$emean
      szpft = datum$szpft
      mdist = datum$lu
      #------------------------------------------------------------------------------------#


      #----- Get soil layers, dslz becomes the weight factor. -----------------------------#
      slz                 = datum$slz
      lsl                 = which.max(x=slz[slz < slz.use])
      slz[lsl]            = slz.use
      dslz                = diff(c(slz,0))
      dslz[slz < slz.use] = 0
      #------------------------------------------------------------------------------------#


      #----- Handy indices. ---------------------------------------------------------------#
      idx   = match(datum$when,when)
      sel   = is.finite(idx)
      #------------------------------------------------------------------------------------#



      #----- Loop over variables. ---------------------------------------------------------#
      for (v in sequence(nvariables)){
         #----- Alias for variable elemnr
         v.vnam   = variables$vnam  [v]
         v.desc   = variables$desc  [v]
         v.unit   = variables$unit  [v]
         v.soil   = variables$soil  [v]
         v.census = variables$census[v]
         v.dset   = variables$dset  [v]
         v.ilu    = variables$ilu   [v]
         v.cvnam  = paste("census",v.vnam,sep=".")
         #---------------------------------------------------------------------------------#


         #----- Grab variables. -----------------------------------------------------------#
         cat ("     # ",v.desc,"...","\n")
         if (v.dset %in% "mdist"){
            areas                = mdist$area[sel,-(nlu+1)]
            lambdas              = mdist$dist[sel,,4]
            output[idx[sel],v,u,s] = log( rowSums(areas) / rowSums(areas*exp(-lambdas)) )
         }else if (v.dset %in% "emean"){
            if (v.soil){
               output[idx[sel],v,u,s] = apply( X      = emean[[v.vnam]][sel,,drop=FALSE]
                                             , MARGIN = 1
                                             , FUN    = weighted.mean
                                             , w      = dslz
                                             )#end apply
            }else{
               output[idx[sel],v,u,s] = emean[[v.vnam]][sel]
            }#end if (v.soil)
         }else{
            if (v.census){
               output[idx[sel] ,v,u,s] = szpft[[v.cvnam]][sel,ndbh+1,npft+1 ]
            }else{
               output[idx[sel] ,v,u,s] = szpft[[v.vnam ]][sel,ndbh+1,npft+1 ]
            }#end if (v.census)
            outsss[idx[sel],,,v,u,s]   = szpft[[v.vnam ]][sel,dbh.use,pft.use]
         }#end if (v.dset %in% "emean")
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(nvariables))
      #------------------------------------------------------------------------------------#
   }#end for (u in sequence(nsimul))
   #---------------------------------------------------------------------------------------#
}#end for (s in sequence(nsites))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Check whether to aggregate data or not.                                              #
#------------------------------------------------------------------------------------------#
if (out.monthly){
   elapsed = whenmonth + 12 * (whenyear - min(whenyear))

   #----- Controls for the X axis. --------------------------------------------------------#
   whenlim = range(elapsed)
   dyr     = mean(diff(pretty(whenlim/12)))
   if (dyr <= 0.1){
      dmon = 1
   }else if (dyr == 0.2){
      dmon = 3
   }else if (dyr == 0.5){
      dmon = 6
   }else{
      dmon = dyr * 12
   }#end if
   whenat = seq(from=0,to=whenlim[2]+dmon,by=dmon)
   whenat = whenat[whenat %wr% whenlim]
   if (dmon >= 12){
      whenlabels = round(whenat / 12)
      el.unit    = untab$yr
   }else{
      whenlabels = whenat
      el.unit    = untab$month
   }#end if
   #---------------------------------------------------------------------------------------#
}else{

   #----- Elapsed time. -------------------------------------------------------------------#
   elapsed = sort(unique(whenyear)) - min(whenyear)
   el.unit = untab$yr
   #---------------------------------------------------------------------------------------#


   #----- Aggregate data by year. ---------------------------------------------------------#
   output  = qapply(X=output,DIM=1,INDEX=whenyear,FUN=mean,na.rm=TRUE)
   outsss  = qapply(X=outsss,DIM=1,INDEX=whenyear,FUN=mean,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#

   #---- Fix dimension names. -------------------------------------------------------------#
   dimnames(output) = list(elapsed,variables$vnam,simul$name,sites$iata)
   dimnames(outsss) = list(elapsed,dbh.key,pft.key,variables$vnam,simul$name,sites$iata)
   #---------------------------------------------------------------------------------------#


   #----- Controls for the X axis. --------------------------------------------------------#
   whenlim    = range(elapsed)
   whenat     = pretty(whenlim)
   whenlabels = whenat
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Plot time series.                                                                    #
#------------------------------------------------------------------------------------------#
for (s in sequence(nsites)){
   s.iata = sites$iata[s]
   s.desc = sites$desc[s]
   cat0(" + ",s.desc," (",toupper(s.iata),")")
   cat0("   - Plot times series by site:")

   #----- Create output path for site. ----------------------------------------------------#
   outsite = file.path(outroot,s.iata)
   if (! file.exists(outsite)) dir.create(outsite)
   #---------------------------------------------------------------------------------------#


   #----- Loop through variables. ---------------------------------------------------------#
   for (v in sequence(nvariables)){
      #----- Variable name. ---------------------------------------------------------------#
      v.vnam = variables$vnam[v]
      v.desc = variables$desc[v]
      v.unit = variables$unit[v]
      v.zero = variables$zero[v]
      cat("   - ",v.desc,"...","\n",sep="")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the range for variables.                                                  #
      #------------------------------------------------------------------------------------#
      if (v.zero){
         vlim    = pretty.xylim(u=c(0,max(output[,v,,s],na.rm=TRUE)))
         vat     = pretty(vlim)
         vlabels = sprintf("%g",vat)
      }else{
         vlim    = pretty.xylim(u=output[,v,,s])
         vat     = pretty(vlim)
         vlabels = sprintf("%g",vat)
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot the time series for all sites.                                            #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Make the file name. -------------------------------------------------------#
         fichier = file.path(outsite,paste0(s.iata,"_tseries_",v.vnam,".",outform[o]))
         if (outform[o] %in% "x11"){
            X11(width=ssize$width,height=ssize$height,pointsize=ptsz)
         }else if (outform[o] %in% "quartz"){
            quartz(width=ssize$width,height=ssize$height,pointsize=ptsz)
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
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Split plot area into two.                                                  #
         #---------------------------------------------------------------------------------#
         par(par.user)
         layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
         #---------------------------------------------------------------------------------#




         #----- Draw legend. --------------------------------------------------------------#
         par(mar=c(0.1,4.1,0.1,1.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x = "bottom"
               , inset = 0.0
               , legend = simul$desc
               , col    = simul$col
               , lty    = simul$lty
               , lwd    = 1.5
               , ncol   = min(3,pretty.box(nsimul)$ncol)
               , cex    = 0.6
               , bty    = "n"
               , xpd    = TRUE
               )#end legend
         #---------------------------------------------------------------------------------#



         #----- Prepare plotting area. ----------------------------------------------------#
         par(mar=c(4.1,4.1,1.6,1.1))
         plot.new()
         plot.window(xlim=whenlim,ylim=vlim)
         axis  (side=1,las=1,at=whenat,labels=whenlabels)
         axis  (side=2,las=1,at=vat   ,labels=vlabels   )
         if (plotgrid) abline(h=yat,v=xat,col=grid.colour,lty="dotted")
         title ( xlab = desc.unit(desc="Simulation Time",unit=el.unit), line = 2.6)
         title ( ylab = desc.unit(desc=v.desc,unit=untab[[v.unit]])    , line = 2.6)
         title ( main = s.desc, line = 0.6, cex.main = 1.0)
         #---------------------------------------------------------------------------------#


         #----- Loop over sites. ----------------------------------------------------------#
         for (u in sequence(nsimul)){
            u.col = simul$col[u]
            u.lty = simul$lty[u]
            lines(x=elapsed,y=output[,v,u,s],col=u.col,lwd=1.5,lty=u.lty)
         }#end for (u in sequence(nsimul))
         #---------------------------------------------------------------------------------#



         #----- Don't forget the box! -----------------------------------------------------#
         box()
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
   }#end for (v in sequence(nvariables))
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Set the x axis for the bar plot sequences.                                        #
   #---------------------------------------------------------------------------------------#
   cat("   - Plot har plot sequences by PFT and DBH:","\n")
   if (out.monthly){
      yplot   = which((elapsed %% (12*ystep)) == 0)
   }else{
      yplot   = which((elapsed %% ystep) == 0)
   }#end if (out.monthly)
   nyplot  = length(yplot)
   xat     = rep(sequence(nyplot)-1,each=ndbh)*(ndbh+1) + rep(sequence(ndbh),times=nyplot)
   xmid    = rep(xat,each=npft.use-1)
   xleft   = xmid - 0.4
   xright  = xmid + 0.4
   xgrid   = (ndbh+1) * (sequence(nyplot+1)-1)
   xlabels = rep( x     = gsub(pattern="cm",replacement="",x=dbh.desc[-ndbh.use])
                , times = nyplot
                )#end rep
   xhead   = mid.points(xgrid)
   xlim    = range(xgrid)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over each variable.                                                          #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(nvariables)){
      #----- Variable name. ---------------------------------------------------------------#
      v.vnam = variables$vnam[v]
      v.desc = variables$desc[v]
      v.unit = variables$unit[v]
      v.dset = variables$dset[v]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Do not plot if this is a "emean" variable.                                     #
      #------------------------------------------------------------------------------------#
      if (v.dset %in% "szpft"){
         #----- Show variable name. -------------------------------------------------------#
         cat("   - ",v.desc,"...","\n",sep="")
         #---------------------------------------------------------------------------------#




         #----- Grab current variable and find the limits for different rectangles. -------#
         outnow  = outsss[yplot,-ndbh.use,-npft.use,v,,s]
         zero    = 0 * outsss[yplot,-ndbh.use,1,v,,s]
         outnow  = abind(zero,outnow,along=3)
         outnow  = apply(X=outnow,MARGIN=c(1,2,4),FUN=cumsum)
         outnow  = aperm(a=outnow,perm=c(1,3,2,4))
         ybottom = outnow[-npft.use,,,]
         ytop    = outnow[       -1,,,]
         dimnames(ybottom)[[1]] = pft.key[-npft.use]
         dimnames(ytop)   [[1]] = pft.key[-npft.use]
         rm(outnow,zero)
         #---------------------------------------------------------------------------------#



         #----- Find the colours. ---------------------------------------------------------#
         idxnow  = arrayInd(seq_along(ytop),.dim=dim(ytop))[,1]
         rectcol = array(data=pft.colour[idxnow],dim=dim(ytop),dimnames=dimnames(ytop))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the Y-axis limits.                                                     #
         #---------------------------------------------------------------------------------#
         ylim    = pretty.xylim(u=c(ytop,ybottom),fracexp=0.1)
         yat     = pretty(ylim)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Loop over simulations.                                                      #
         #---------------------------------------------------------------------------------#
         for (u in sequence(nsimul)){
            u.name = simul$name[u]
            u.desc = simul$desc[u]
            cat("     * ",u.desc,"...","\n",sep="")



            #------------------------------------------------------------------------------#
            #     Plot the snapshots of each variable.                                     #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Make the file name. -------------------------------------------------#
               fichier = file.path(outsite
                                  ,paste0(s.iata,"_",u.name,"_barplot-",v.vnam,"-"
                                         ,".",outform[o]))
               if (outform[o] %in% "x11"){
                  X11(width=wsize$width,height=wsize$height,pointsize=ptsz)
               }else if (outform[o] %in% "quartz"){
                  quartz(width=wsize$width,height=wsize$height,pointsize=ptsz)
               }else if(outform[o] %in% "png"){
                  png(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
                     ,pointsize=ptsz,res=depth,bg="transparent")
               }else if(outform[o] %in% "tif"){
                  tiff(filename=fichier,width=wsize$width*depth,height=wsize$height*depth
                      ,pointsize=ptsz,res=depth,bg="transparent",compression="lzw")
               }else if(outform[o] %in% "eps"){
                  postscript(file=fichier,width=wsize$width,height=wsize$height
                            ,pointsize=ptsz,paper=wsize$paper)
               }else if(outform[o] %in% "pdf"){
                  pdf(file=fichier,onefile=FALSE,width=wsize$width,height=wsize$height
                     ,pointsize=ptsz,paper=wsize$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split device area into two.                                           #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     First plot, the legend.                                               #
               #---------------------------------------------------------------------------#
               par(mar=c(0.1,4.1,0.1,1.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x      = "bottom"
                     , inset  = 0.0
                     , legend = pft.desc  [-npft.use]
                     , border = pft.colour[-npft.use]
                     , fill   = pft.colour[-npft.use]
                     , horiz  = TRUE
                     , title  = expression(bold("Plant functional types"))
                     , cex    = 0.8
                     , bty    = "n"
                     , xpd    = TRUE
                     )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Plot the bars.                                                        #
               #---------------------------------------------------------------------------#
               par(mar=c(4.1,4.1,3.1,1.1))
               plot.new()
               plot.window(xlim=xlim,ylim=ylim)
               axis(side=1,las=2,at=xat,labels=xlabels,cex.axis=0.8)
               axis(side=2,las=1,at=yat,labels=ylabels)
               title(xlab= desc.unit(desc="DBH class",unit=untab$cm)  , cex.lab  = 0.8)
               title(ylab= desc.unit(desc=v.desc,unit=untab[[v.unit]]), line     = 2.6)
               title(main= paste0(s.desc," - ",u.desc) , cex.main = 1.0)
               abline(v=xgrid,col=grey.fg,lwd=0.8,lty="dotdash")
               abline(h=0    ,col=grey.fg,lwd=0.8,lty="solid"  )
               rect( xleft   = xleft
                   , xright  = xright
                   , ybottom = ybottom[,,,u]
                   , ytop    = ytop   [,,,u]
                   , col     = rectcol[,,,u]
                   , density = NA
                   , lwd     = 2
                   )#end rect
               text( x      = xhead
                   , y      = ylim[1]+0.98*diff(ylim)
                   , labels = paste0("+",elapsed[yplot])
                   , adj    = c(0.5,1)
                   )#end text
               box()
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
            }#end for (o in sequence(nout)){
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsites))
         #---------------------------------------------------------------------------------#
      }#end if (v.dset %in% "szpft")
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(nvariables))
   #---------------------------------------------------------------------------------------#
}#end for (s in sequence(nsites))
#------------------------------------------------------------------------------------------#


