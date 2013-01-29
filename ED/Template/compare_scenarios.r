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
#------------------------------------------------------------------------------------------#
here    = getwd()                                  #   Current directory
srcdir  = "/n/moorcroft_data/mlongo/util/Rsc"      #   Script directory
outroot = paste(here,"scenario_evergreen",sep="/") #   Output directory

sites      = c("gyf","s67","m34","rja","cax","pdg","ban","pnz")
simul      = list()
simul[[1]] = list( name   = "iscen-01_phen-01_canrad01"
                 , key    = "past"
                 , desc   = "Historical rainfall"
                 , colour = "royalblue4"
                 , pch    = 15
                 )#end list
simul[[2]] = list( name   = "iscen+01_phen-01_canrad01"
                 , key    = "unif"
                 , desc   = "Uniform resampling"
                 , colour = "steelblue3"
                 , pch    = 17
                 )#end list
simul[[3]] = list( name   = "iscen-02_phen-01_canrad01"
                 , key    = "dry02"
                 , desc   = "Dry bias (A = -2)"
                 , colour = "orange1"
                 , pch    =  9
                 )#end list
simul[[4]] = list( name   = "iscen-06_phen-01_canrad01"
                 , key    = "dry06"
                 , desc   = "Dry bias (A = -6)"
                 , colour = "firebrick"
                 , pch    = 13
                 )#end list
#------------------------------------------------------------------------------------------#





#------ Miscellaneous settings. -----------------------------------------------------------#
yeara          = 1972         # First year we will include
yearz          = 2011         # Last year we will include
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      = 2            # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
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
ptsz           = 14                    # Font size.
lwidth         = 2.5                   # Line width
plotgrid       = TRUE                  # Should I plot the grid in the background? 

legwhere       = "topleft"             # Where should I place the legend?
inset          = 0.01                  # Inset between legend and edge of plot region.
legbg          = "white"               # Legend background colour.
fracexp        = 0.40                  # Expansion factor for y axis (to fit legend)
cex.main       = 0.8                   # Scale coefficient for the title
xyz.ncolour    = 20                    # Number of colours for the xyz
notch          = FALSE                 # Add notches to the box plots.
mtext.xoff     = -7.00                 # Offset for the x label
mtext.yoff     = -1.00                 # Offset for the y label
mtext.xadj     =  0.50                 # Offset for the x label
mtext.yadj     =  0.65                 # Offset for the y label
barplot.lwd    =  2.00                 # Line width for the bar plots
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



#----- Set some dimensions associated with the simulations. -------------------------------#
n.sites = length(sites)
n.simul = length(simul)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Replace the list by a data frame.                                                    #
#------------------------------------------------------------------------------------------#
simul = data.frame( apply( X = sapply(X=simul,FUN=c), MARGIN = 1, FUN = unlist )
                  , stringsAsFactors = FALSE
                  )#end data.frame
for (nn in c("pch")) simul[[nn]] = as.numeric(simul[[nn]])
#------------------------------------------------------------------------------------------#




#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#      Create the output paths.                                                            #
#------------------------------------------------------------------------------------------#

   #----- Make sure that the base directory exists. ---------------------------------------#
   if (! file.exists(outroot)) dir.create(outroot)
   #---------------------------------------------------------------------------------------#




   #----- Create the paths by type of plot. -----------------------------------------------#
   root.ts.season      = file.path(outroot,"ts_season"     )
   root.ts.year        = file.path(outroot,"ts_year"       )
   root.tspft.season   = file.path(outroot,"tspft_season"  )
   root.tspft.year     = file.path(outroot,"tspft_year"    )
   root.boxpft.season  = file.path(outroot,"boxpft_season" )
   root.boxpft.year    = file.path(outroot,"boxpft_year"   )
   root.boxpftdbh      = file.path(outroot,"boxpftdbh"     )
   root.barplot.season = file.path(outroot,"barplot_season")
   root.barplot.year   = file.path(outroot,"barplot_year"  )
   root.xyz.season     = file.path(outroot,"xyz_season"    )
   root.xyz.pft        = file.path(outroot,"xyz_pft"       )
   if (! file.exists(root.ts.season     )) dir.create(root.ts.season     )
   if (! file.exists(root.ts.year       )) dir.create(root.ts.year       )
   if (! file.exists(root.tspft.season  )) dir.create(root.tspft.season  )
   if (! file.exists(root.tspft.year    )) dir.create(root.tspft.year    )
   if (! file.exists(root.boxpft.season )) dir.create(root.boxpft.season )
   if (! file.exists(root.boxpft.year   )) dir.create(root.boxpft.year   )
   if (! file.exists(root.boxpftdbh     )) dir.create(root.boxpftdbh     )
   if (! file.exists(root.barplot.season)) dir.create(root.barplot.season)
   if (! file.exists(root.barplot.year  )) dir.create(root.barplot.year  )
   if (! file.exists(root.xyz.season    )) dir.create(root.xyz.season    )
   if (! file.exists(root.xyz.pft       )) dir.create(root.xyz.pft       )
   #---------------------------------------------------------------------------------------#




   #----- Generate the names of the sub-sub-directories by variable. ----------------------#
   out.ts.season      = file.path(root.ts.season     ,scen.ts$vname      )
   out.ts.year        = file.path(root.ts.year       ,scen.ts$vname      )
   out.tspft.season   = file.path(root.tspft.season  ,scen.ts$vname      )
   out.tspft.year     = file.path(root.tspft.year    ,scen.ts$vname      )
   out.boxpft.season  = file.path(root.boxpft.season ,scen.szpft$vname   )
   out.boxpft.year    = file.path(root.boxpft.year   ,scen.szpft$vname   )
   out.boxpftdbh      = file.path(root.boxpftdbh     ,scen.szpft$vname   )
   out.barplot.season = file.path(root.barplot.season,scen.barplot$vname )
   out.barplot.year   = file.path(root.barplot.year  ,scen.barplot$vname )
   out.xyz.season     = file.path(root.xyz.season    ,scen.xyz$yvar$vname)
   out.xyz.pft        = file.path(root.xyz.pft       ,scen.xyz$yvar$vname)
   #---------------------------------------------------------------------------------------#



   #----- Loop over all time series variables and create the paths as needed. -------------#
   for (n in 1:nscen.ts){
      #----- Check whether this variable is to be plotted. --------------------------------#
      if (scen.ts$plt[n]){
         if (! file.exists(out.ts.season[n])) dir.create(out.ts.season[n])
         if (! file.exists(out.ts.year  [n])) dir.create(out.ts.year  [n])
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Check whether this variable is to be plotted by PFT. -------------------------#
      if (scen.ts$plt[n] && scen.ts$pftvar[n]){
         if (! file.exists(out.tspft.season[n])) dir.create(out.tspft.season[n])
         if (! file.exists(out.tspft.year  [n])) dir.create(out.tspft.year  [n])
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Loop over all box plot variables and create the paths as needed. ----------------#
   for (n in 1:nscen.szpft){
      #----- Check whether this variable is to be plotted. --------------------------------#
      if (scen.szpft$plt[n]){
         if (! file.exists(out.boxpft.season[n])) dir.create(out.boxpft.season[n])
         if (! file.exists(out.boxpft.year  [n])) dir.create(out.boxpft.year  [n])
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Check whether this variable is to be plotted by PFT. -------------------------#
      if (scen.szpft$plt[n] && scen.szpft$dbhvar[n]){
         if (! file.exists(out.boxpftdbh[n])) dir.create(out.boxpftdbh[n])
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Loop over all bar plot variables and create the paths as needed. ----------------#
   for (n in 1:nscen.barplot){
      #----- Check whether this variable is to be plotted. --------------------------------#
      if (scen.barplot$plt[n]){
         if (! file.exists(out.barplot.season[n])) dir.create(out.barplot.season[n])
         if (! file.exists(out.barplot.year  [n])) dir.create(out.barplot.year  [n])
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Loop over all parameter space variables (Y) and create the paths as needed. -----#
   for (n in 1:nscen.yvar){
      if (! file.exists(out.xyz.season[n])) dir.create(out.xyz.season[n])
      if (! file.exists(out.xyz.pft   [n])) dir.create(out.xyz.pft   [n])
   }#end for
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#



#----- Load observations. -----------------------------------------------------------------#
obsrfile = paste(srcdir,"LBA_MIP.v8.RData",sep="/")
load(file=obsrfile)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define some dimensions.                                                              #
#------------------------------------------------------------------------------------------#
pft.use       = c( 1, 2, 3, 4,16, 18)    # PFT classes to include (add PFT=18, the total)
pft.mp        = c( F, T, T, T, F,  T)    # Include the PFT on multi-panel plots?
pft.key       = pft$key    [pft.use   ]  # PFT keys  (for dimnames)
pft.desc      = pft$name   [pft.use   ]  # PFT names (for titles)
pft.colour    = pft$colour [pft.use   ]  # PFT colours
dbh.use       = seq(from=2,to=ndbh,by=1) # DBH classes that we will use
dbh.key       = dbhkeys    [dbh.use   ]  # DBH keys  (for dimnames)
dbh.desc      = dbhnames   [dbh.use   ]  # DBH names (for titles)
dbh.colour    = dbhcols    [dbh.use   ]  # DBH colours
season.use    = c(1,2,3,4,5)             # Seasons to include (add season=5, the total)
season.key    = season.list[season.use]  # Keys for the seasons. 
season.desc   = season.full[season.use]  # Full names of all seasons.
season.colour = season.cols[season.use]  # Colours for seasons.
year.use      = yeara:yearz              # Years to use
year.key      = year.use                 # Year keys  (for dimnames)
year.desc     = year.use                 # Year names (for titles)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Transform pft.mp into an index linked to pft.use and associated variables.           #
#------------------------------------------------------------------------------------------#
pft.mp           = which(pft.mp)
pft.bp           = sequence(length(pft.use)-1)
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
n.season      = length(season.use)       # Number of seasons
n.season.mp   = n.season-1               # Number of seasons for multiple panels
n.year        = length(year.use)         # Number of years
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Find the best set up for plotting all seasons and all PFTs in the same plot.        #
#------------------------------------------------------------------------------------------#
lo.season = pretty.box(n=n.season.mp)
lo.pft    = pretty.box(n=n.pft.mp   )
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Loop over all sites.                                                                 #
#------------------------------------------------------------------------------------------#
for (p in 1:n.sites){
   #----- Retrieve the site information. --------------------------------------------------#
   pidx     = match(sites[p],poilist$iata)
   short    = poilist$short   [pidx]
   longname = poilist$longname[pidx]
   iata     = poilist$iata    [pidx]
   lon      = poilist$lon     [pidx]
   lat      = poilist$lat     [pidx]
   cat(" + Comparing simulations for ",longname,"...","\n",sep="")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Initialise the list of variables.                                                 #
   #---------------------------------------------------------------------------------------#
   eft = list()
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create a general time stamp that works for all simulations.                       #
   #---------------------------------------------------------------------------------------#
   #----- Get all times. ------------------------------------------------------------------#
   eft$when      = c(chron(paste(12,1,yeara-1,sep="/"))
                    ,chron(paste(rep(1:12,times=n.year),1,rep(year.use,each=12),sep="/"))
                    )#end c
   eft$when      = eft$when[-length(eft$when)]
   eft$year      = numyears (eft$when)
   eft$month     = nummonths(eft$when)
   n.when        = length(eft$when)
   #---- Find the seasons. ----------------------------------------------------------------#
   eft$season    = season(eft$when,add.year=TRUE,dec.next=TRUE)
   eft$ss.year   = as.numeric(substring(eft$season,1,4))
   eft$ss.season = as.numeric(substring(eft$season,5,6))
   #---- Find unique identifiers for seasons and years based on season, not months. -------#
   eft$toseason  = unique(eft$season)
   eft$toyear    = unique(eft$ss.year)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make some arrays with the right dimensions.                                       #
   #---------------------------------------------------------------------------------------#
   empty          = rep  (NA,times=n.when)
   empty.pft      = array(NA,dim=c(n.when,n.pft))
   empty.pftdbh   = array(NA,dim=c(n.when,n.dbh,n.pft))
   ts.array       = array( data     = NA
                         , dim      = c   (  n.year,  n.season,  n.simul)
                         , dimnames = list(year.key,season.key,simul$key)
                         )#end array
   tspft.array    = array( data     = NA
                         , dim      = c   (  n.year,  n.season,  n.pft,  n.simul)
                         , dimnames = list(year.key,season.key,pft.key,simul$key)
                         )#end array
   tspftdbh.array = array( data     = NA
                         , dim      = c   (   n.year,  n.season,  n.dbh,  n.pft,  n.simul)
                         , dimnames = list( year.key,season.key,dbh.key,pft.key,simul$key)
                         )#end array
   #---------------------------------------------------------------------------------------#




   #=======================================================================================#
   #=======================================================================================#
   #     Initialise all variables.                                                         #
   #---------------------------------------------------------------------------------------#
   for (v in 1:nscen.ts){
      #----- Copy variable info. ----------------------------------------------------------#
      var.vname  = scen.ts$vname[v]
      var.desc   = scen.ts$desc [v]
      var.quant  = scen.ts$quant[v]
      is.pft     = scen.ts$pft  [v]
      is.dbh     = scen.ts$dbh  [v]
      var.pft    = paste(var.vname,"pft"   ,sep="")
      var.pftdbh = paste(var.vname,"pftdbh",sep="")
      cat  ("     * Creating data holders for ",var.desc,"...","\n")

      #----- Append the lists to a common name. -------------------------------------------#
      eft[[var.vname]]                      = list()
      eft[[var.vname]]$ts                   = ts.array
      if (is.pft) eft[[var.vname]]$tspft    = tspft.array
      if (is.dbh) eft[[var.vname]]$tspftdbh = tspftdbh.array
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:nscen.ts)
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     Loop over the simulations.                                                        #
   #---------------------------------------------------------------------------------------#
   for (s in 1:n.simul){
      #---- Copy the settings to local variables. -----------------------------------------#
      sim.name   = simul$name  [s]
      sim.desc   = simul$desc  [s]
      sim.colour = simul$colour[s]
      sim.pch    = simul$pch   [s]
      sim.full   = paste("t",iata,"_",sim.name,sep="")
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Load the data set.                                                             #
      #------------------------------------------------------------------------------------#
      rdata.simul = paste(here,"/",sim.full,"/rdata_month/",sim.full,".RData",sep="")
      cat  ("   - Load data from file ",basename(rdata.simul),"...","\n")
      load (rdata.simul) 
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the indices for mapping the data from the original data set to the        #
      # combined one.                                                                      #
      #------------------------------------------------------------------------------------#
      idx = match(eft$when,datum$when)
      sel = is.finite(idx)
      idx = idx[sel]
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      These are shorter versions of the season indices.                             #
      #------------------------------------------------------------------------------------#
      ee = sequence(n.season.mp)
      e5 = n.season
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Copy the time variables to the consolidated list.                              #
      #------------------------------------------------------------------------------------#
      for (v in 1:nscen.ts){
         #----- Copy variable info. -------------------------------------------------------#
         var.vname  = scen.ts$vname [v]
         var.desc   = scen.ts$desc  [v]
         var.f.aggr = get(scen.ts$f.aggr[v])
         is.pft     = scen.ts$pft   [v]
         is.dbh     = scen.ts$dbh   [v]
         is.mort    = scen.ts$mort  [v]
         is.recr    = scen.ts$recr  [v]
         var.pft    = paste(var.vname,"pft"   ,sep="")
         var.pftdbh = paste(var.vname,"pftdbh",sep="")
         cat  ("     * Processing ",var.desc,"...","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Grab the time series.                                                       #
         #---------------------------------------------------------------------------------#
         var.now    = empty
         if (is.pft){
            var.now[sel] = datum[[var.pft]][idx,pft.use[n.pft]]
         }else{
            var.now[sel] = datum[[var.vname]][idx]
         }#end if (is.pft)
         #------ Find the means/sum by year and by season. --------------------------------#
         eft[[var.vname]]$ts[,ee,s] = tapply( X     = var.now
                                            , INDEX = list(eft$ss.year,eft$ss.season)
                                            , FUN   = var.f.aggr
                                            , na.rm = TRUE
                                            )#end tapply
         eft[[var.vname]]$ts[,e5,s] = tapply( X     = var.now
                                            , INDEX = eft$ss.year
                                            , FUN   = var.f.aggr
                                            , na.rm = TRUE
                                            )#end tapply
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check what to do depending on whether the variable is a PFT and/or DBH.     #
         #---------------------------------------------------------------------------------#
         if (is.pft){
            var.now       = empty.pft
            var.now[sel,] = datum[[var.pft]][idx,pft.use]

            #----- Find the means/sum by year and by season. ------------------------------#
            eft[[var.vname]]$tspft[,ee,,s] = qapply( X     = var.now
                                                   , INDEX = list( eft$ss.year
                                                                 , eft$ss.season
                                                                 )#end list
                                                   , DIM   = 1
                                                   , FUN   = var.f.aggr
                                                   , na.rm = TRUE
                                                   )#end qapply
            eft[[var.vname]]$tspft[,e5,,s] = qapply( X     = var.now
                                                   , INDEX = eft$ss.year
                                                   , DIM   = 1
                                                   , FUN   = var.f.aggr
                                                   , na.rm = TRUE
                                                   )#end qapply
            #------------------------------------------------------------------------------#
         }#end if (var.pft)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check what to do depending on whether the variable is a PFT and/or DBH.     #
         #---------------------------------------------------------------------------------#
         if (is.dbh){
            var.now        = empty.pftdbh
            var.now[sel,,] = datum[[var.pftdbh]][idx,dbh.use,pft.use]

            #----- Find the means/sum by year and by season. ------------------------------#
            eft[[var.vname]]$tspftdbh[,ee,,,s] = qapply( X     = var.now
                                                       , INDEX = list( eft$ss.year
                                                                     , eft$ss.season
                                                                     )#end list
                                                       , DIM   = 1
                                                       , FUN   = var.f.aggr
                                                       , na.rm = TRUE
                                                       )#end qapply
            eft[[var.vname]]$tspftdbh[,e5,,,s] = qapply( X     = var.now
                                                       , INDEX = eft$ss.year
                                                       , DIM   = 1
                                                       , FUN   = var.f.aggr
                                                       , na.rm = TRUE
                                                       )#end qapply
            #------------------------------------------------------------------------------#
         }#end if (var.dbh)
         #---------------------------------------------------------------------------------#
      }#end for (v in 1:nscen.ts)
      #------------------------------------------------------------------------------------#
   }#end for (s in 1:n.simul)
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      Loop over variables and find the ones that are mortality or recruitment.  We     #
   # must transform them so they are always between 0 and 100.                             #
   #---------------------------------------------------------------------------------------#
   cat  ("     * Transforming mortality and recruitment variables...","\n")
   for (v in 1:nscen.ts){
      #----- Copy variable info. ----------------------------------------------------------#
      var.vname  = scen.ts$vname[v]
      var.desc   = scen.ts$desc [v]
      is.pft     = scen.ts$pft  [v]
      is.dbh     = scen.ts$dbh  [v]
      is.mort    = scen.ts$mort [v]
      is.recr    = scen.ts$recr [v]
      cat  ("     * Creating data holders for ",var.desc,"...","\n")

      #----- Transform mortality data. ----------------------------------------------------#
      if (is.mort){
         eft[[var.vname]]$ts                  = 100.*(1. - exp(-eft[[var.vname]]$ts      ))
         if(is.pft) eft[[var.vname]]$tspft    = 100.*(1. - exp(-eft[[var.vname]]$tspft   ))
         if(is.dbh) eft[[var.vname]]$tspftdbh = 100.*(1. - exp(-eft[[var.vname]]$tspftdbh))
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Transform recruitment data. --------------------------------------------------#
      if (is.recr){
         eft[[var.vname]]$ts                  = 100.*(exp(eft[[var.vname]]$ts      ) - 1.0)
         if(is.pft) eft[[var.vname]]$tspft    = 100.*(exp(eft[[var.vname]]$tspft   ) - 1.0)
         if(is.dbh) eft[[var.vname]]$tspftdbh = 100.*(exp(eft[[var.vname]]$tspftdbh) - 1.0)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:nscen.ts)
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     Plot the comparison between time series, by season.                               #
   #---------------------------------------------------------------------------------------#
   cat  ("   - Plotting the time series...","\n")
   for (v in 1:nscen.ts){
      #----- Copy variable info. ----------------------------------------------------------#
      var.vname  = scen.ts$vname [v]
      var.desc   = scen.ts$desc  [v]
      var.unit   = scen.ts$unit  [v]
      is.pft     = scen.ts$pftvar[v]
      var.plog   = scen.ts$plog  [v]
      var.plt    = scen.ts$plt   [v]
      if (var.plog){
         var.xylog  = "y"
      }else{
         var.xylog  = ""
      }#end if

      if (var.plt){
         cat  ("   - ",var.desc,"...","\n")
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the variable by season.                                               #
         #---------------------------------------------------------------------------------#


         #----- Load variable. ------------------------------------------------------------#
         this    = eft[[var.vname]]$ts
         xlimit  = pretty.xylim(year.use             ,fracexp=0.0,is.log=FALSE   )
         ylimit  = pretty.xylim(this[,1:n.season.mp,],fracexp=0.0,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Seasonal means) - ",longname,sep="")
         lex     = paste("Year")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.ts.season[v],"/",var.vname,"-",iata,"-ts_season."
                           ,outform[o],sep="")
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
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,4,0))
            layout(mat    = rbind(1+lo.season$mat,rep(1,times=lo.season$ncol))
                  ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , col     = simul$colour
                   , lwd     = 2.0
                   , pch     = 16
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all seasons, and plot the bar plots.                           #
            #------------------------------------------------------------------------------#
            for (e in 1:n.season.mp){
               #----- Find out where is this box going, and set up axes and margins. ------#
               left    = (e %% lo.season$ncol) == 1
               right   = (e %% lo.season$ncol) == 0
               top     = e <= lo.season$ncol
               bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
               mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
               #---------------------------------------------------------------------------#


               #----- Set up the title for each plot. -------------------------------------#
               lesub = paste(season.desc[e],sep="")
               #---------------------------------------------------------------------------#



               #----- Plot window and grid. -----------------------------------------------#
               par(mar=mar.now)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
               if (bottom) axis(side=1)
               if (left  ) axis(side=2)
               box()
               title(main=lesub)
               if (plotgrid) grid(col="grey83",lty="solid")

               #----- Plot the lines. -----------------------------------------------------#
               for (s in 1:n.simul){
                  points(x=year.use,y=this[,e,s],type="o",pch=16,col=simul$colour[s]
                        ,lwd=2.0)
               }#end for
               #---------------------------------------------------------------------------#
            }#end for (e in 1:n.season.mp)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
            mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Now we plot the annual data.                                                #
         #---------------------------------------------------------------------------------#


         #----- Update limits. ------------------------------------------------------------#
         xlimit  = pretty.xylim(year.use        ,fracexp=0.0    ,is.log=FALSE   )
         ylimit  = pretty.xylim(this[,n.season,],fracexp=fracexp,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Annual means) ",longname,sep="")
         lex     = paste("Year")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.ts.year[v],"/",var.vname,"-",iata,"-ts_year."
                           ,outform[o],sep="")
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



            #----- Plot window and grid. --------------------------------------------------#
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            axis(side=1)
            axis(side=2)
            box()
            title(main=letitre,xlab=lex,ylab=ley)
            if (plotgrid) grid(col="grey83",lty="solid")
            #------------------------------------------------------------------------------#



            #----- Plot the lines. --------------------------------------------------------#
            for (s in 1:n.simul){
               points(x=year.use,y=this[,n.season,s],type="o",pch=16,col=simul$colour[s]
                     ,lwd=2.0)
            }#end for
            #------------------------------------------------------------------------------#



            #---- Add the legend. ---------------------------------------------------------#
            legend ( x       = "topright"
                   , inset   = inset
                   , legend  = simul$desc
                   , col     = simul$colour
                   , lwd     = 2.0
                   , pch     = 16
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end if (var.plt)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     If this is a PFT variable, plot the time series by PFT and season.             #
      #------------------------------------------------------------------------------------#
      if (var.plt && is.pft){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #     Plot the seasonal time series for each PFT.                                 #
         #---------------------------------------------------------------------------------#
         for (f in 1:(n.pft-1)){
            cpft = sprintf("%2.2i",pft.use[f])

            #----- Load variable. ---------------------------------------------------------#
            this    = eft[[var.vname]]$tspft
            xlimit  = pretty.xylim(year.use               ,fracexp=0.0,is.log=FALSE   )
            ylimit  = pretty.xylim(this[,1:n.season.mp,f,],fracexp=0.0,is.log=var.plog)
            #------------------------------------------------------------------------------#



            #----- Set the title. ---------------------------------------------------------#
            letitre = paste(var.desc," (Seasonal means) - ",longname
                           ,"\n",pft.desc[f],sep="")
            lex     = paste("Year")
            ley     = paste(var.desc," [",var.unit,"]",sep=" ")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all formats.                                                  #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open file or display. -----------------------------------------------#
               fichier = paste(out.tspft.season[v],"/",var.vname,"-",iata,"-",cpft
                              ,"-ts_season.",outform[o],sep="")
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
               #     Split the window into several smaller windows.  Add a bottom row to   #
               # fit the legend.                                                           #
               #---------------------------------------------------------------------------#
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,4,0))
               layout(mat    = rbind(1+lo.season$mat,rep(1,times=lo.season$ncol))
                     ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#



               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = simul$desc
                      , col     = simul$colour
                      , lwd     = 2.0
                      , pch     = 16
                      , bg      = "white"
                      , ncol    = min(3,pretty.box(n.simul)$ncol)
                      , title   = expression(bold("Simulation"))
                      , cex     = 1.0
                      )#end legend
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over all seasons, and plot the bar plots.                        #
               #---------------------------------------------------------------------------#
               for (e in 1:n.season.mp){
                  #----- Find out where is this box going, and set up axes and margins. ---#
                  left    = (e %% lo.season$ncol) == 1
                  right   = (e %% lo.season$ncol) == 0
                  top     = e <= lo.season$ncol
                  bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
                  mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
                  #------------------------------------------------------------------------#


                  #----- Set up the title for each plot. ----------------------------------#
                  lesub = paste(season.desc[e],sep="")
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=mar.now)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  if (bottom) axis(side=1)
                  if (left  ) axis(side=2)
                  box()
                  title(main=lesub)
                  if (plotgrid) grid(col="grey83",lty="solid")

                  #----- Plot the lines. --------------------------------------------------#
                  for (s in 1:n.simul){
                     points(x=year.use,y=this[,e,f,s],type="o",pch=16,col=simul$colour[s]
                           ,lwd=2.0)
                  }#end for
                  #------------------------------------------------------------------------#
               }#end for (e in 1:n.season.mp)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
               mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
               mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end for (f in 1:(n.pft-1))
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the time series by PFT (annual, one panel for each PFT).              #
         #---------------------------------------------------------------------------------#


         #----- Reset limits. -------------------------------------------------------------#
         this    = eft[[var.vname]]$tspft
         xlimit  = pretty.xylim(year.use               ,fracexp=0.0,is.log=FALSE   )
         ylimit  = pretty.xylim(this[,n.season,pft.mp,],fracexp=0.0,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Annual means) - ",longname,sep="")
         lex     = paste("Year")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.tspft.year[v],"/",var.vname,"-",iata,"-tspft_year."
                           ,outform[o],sep="")
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
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,4,0))
            layout(mat    = rbind(1+lo.pft$mat,rep(1,times=lo.pft$ncol))
                  ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , col     = simul$colour
                   , lwd     = 2.0
                   , pch     = 16
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all seasons, and plot the bar plots.                           #
            #------------------------------------------------------------------------------#
            for (f in 1:n.pft.mp){

               #----- Find out where is this box going, and set up axes and margins. ------#
               left    = (f %% lo.pft$ncol) == 1
               right   = (f %% lo.pft$ncol) == 0
               top     = f <= lo.pft$ncol
               bottom  = f > (lo.pft$nrow - 1) * lo.pft$ncol
               mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
               #---------------------------------------------------------------------------#


               #----- Set up the title for each plot. -------------------------------------#
               lesub = paste(pft.desc[pft.mp[f]],sep="")
               #---------------------------------------------------------------------------#



               #----- Plot window and grid. -----------------------------------------------#
               par(mar=mar.now)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
               if (bottom) axis(side=1)
               if (left  ) axis(side=2)
               box()
               title(main=lesub)
               if (plotgrid) grid(col="grey83",lty="solid")

               #----- Plot the lines. -----------------------------------------------------#
               for (s in 1:n.simul){
                  points(x=year.use,y=this[,n.season,pft.mp[f],s],type="o",pch=16
                        ,col=simul$colour[s],lwd=2.0)
               }#end for
               #---------------------------------------------------------------------------#
            }#end for (e in 1:n.pft)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
            mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end if (is.pft && var.plt)
      #------------------------------------------------------------------------------------#

   }#end for (v in 1:nscen.ts)
   #=======================================================================================#
   #=======================================================================================#








   #=======================================================================================#
   #=======================================================================================#
   #      Plot the box plots.                                                              #
   #---------------------------------------------------------------------------------------#
   cat  ("   - Plotting the box plots...","\n")
   for (v in 1:nscen.szpft){
      #----- Copy variable info. ----------------------------------------------------------#
      var.vname  = scen.szpft$vname [v]
      var.desc   = scen.szpft$desc  [v]
      var.unit   = scen.szpft$unit  [v]
      is.pft     = scen.szpft$pftvar[v]
      is.dbh     = scen.szpft$dbhvar[v]
      var.plog   = scen.szpft$plog  [v]
      var.plt    = scen.szpft$plt   [v]
      if (var.plog){
         var.xylog  = "y"
      }else{
         var.xylog  = ""
      }#end if
      #------------------------------------------------------------------------------------#


      if (var.plt){
         cat  ("   - ",var.desc,"...","\n")
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the box plot by PFT and by season.                                    #
         #---------------------------------------------------------------------------------#
         xlimit = c(0,n.pft.bp*n.simul) + 0.5
         xat    = seq(from=1 + 0.5*(n.simul-1), to=n.pft.bp*n.simul, by=n.simul)
         xgrid  = seq(from=0.5, to=(n.pft.bp+1)*n.simul - 0.5*(n.simul-1), by=n.simul)
         #---------------------------------------------------------------------------------#


         #----- Load variable and set y axis. ---------------------------------------------#
         this      = eft[[var.vname]]$tspft[,1:n.season.mp,pft.bp,]
         ylimit    = pretty.xylim(u=this,fracexp=0.0,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Seasonal means) - ",longname,sep="")
         lex     = paste("Plant functional type")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Set up the box plot colour.                                                 #
         #---------------------------------------------------------------------------------#
         bp.colour = rep(simul$colour,times=n.pft.bp)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.boxpft.season[v],"/",var.vname,"-",iata,"-boxplot_season."
                           ,outform[o],sep="")
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
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,4,0))
            layout(mat    = rbind(1+lo.season$mat,rep(1,times=lo.season$ncol))
                  ,height = c(rep(5/lo.season$nrow,times=lo.season$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#



            #----- Plot legend. -----------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , fill    = simul$colour
                   , border  = "black"
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Loop over all seasons, and plot the bar plots.                           #
            #------------------------------------------------------------------------------#
            for (e in 1:n.season.mp){

               #----- Find out where is this box going, and set up axes and margins. ------#
               left    = (e %% lo.season$ncol) == 1
               right   = (e %% lo.season$ncol) == 0
               top     = e <= lo.season$ncol
               bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
               mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
               #---------------------------------------------------------------------------#


               #----- Set up the title for each plot. -------------------------------------#
               lesub = paste(season.desc[e],sep="")
               #---------------------------------------------------------------------------#



               #----- Create the temporary box plot list. ---------------------------------#
               ss.this = this[,e,,]
               ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
               bp.plot = split(x=ss.this,f=list(ai[,3],ai[,2]))
               #---------------------------------------------------------------------------#



               #----- Plot window and grid. -----------------------------------------------#
               par(mar=mar.now)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
               if (bottom) axis(side=1,at=xat,labels=pft.key[pft.bp])
               if (left  ) axis(side=2)
               box()
               title(main=lesub)
               if (plotgrid) abline(v=xgrid,h=axTicks(2),col="grey40",lty="solid")

               #----- Add the box plot, without the x axis. -------------------------------#
               boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                      ,yaxt="n")
               #---------------------------------------------------------------------------#
            }#end for (e in 1:n.season.mp)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
            mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the box plot by PFT and by season.                                    #
         #---------------------------------------------------------------------------------#
         xlimit = c(0,n.pft.bp*n.simul) + 0.5
         xat    = seq(from=1 + 0.5*(n.simul-1), to=n.pft.bp*n.simul, by=n.simul)
         xgrid  = seq(from=0.5, to=(n.pft.bp+1)*n.simul - 0.5*(n.simul-1), by=n.simul)
         #---------------------------------------------------------------------------------#


         #----- Load variable and set y axis. ---------------------------------------------#
         this      = eft[[var.vname]]$tspft[,n.season,pft.bp,]
         ylimit    = pretty.xylim(u=this,fracexp=fracexp,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Annual means) - ",longname,sep="")
         lex     = paste("Plant functional type")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Set up the box plot colour.                                                 #
         #---------------------------------------------------------------------------------#
         bp.colour = rep(simul$colour,times=n.pft.bp)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.boxpft.year[v],"/",var.vname,"-",iata,"-boxplot_year."
                           ,outform[o],sep="")
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



            #----- Create the temporary box plot list. ------------------------------------#
            ai      = arrayInd(sequence(length(this)),.dim=dim(this))
            bp.plot = split(x=ss.this,f=list(ai[,3],ai[,2]))
            #------------------------------------------------------------------------------#



            #----- Plot window and grid. --------------------------------------------------#
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
            if (bottom) axis(side=1,at=xat,labels=pft.key[pft.bp])
            if (left  ) axis(side=2)
            box()
            title(main=lesub)
            if (plotgrid) abline(v=xgrid,h=axTicks(2),col="grey40",lty="solid")

            #----- Add the box plot, without the x axis. ----------------------------------#
            boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE,yaxt="n")
            #------------------------------------------------------------------------------#



            #----- Add the legend. --------------------------------------------------------#
            legend ( x       = "topright"
                   , inset   = inset
                   , legend  = simul$desc
                   , fill    = simul$colour
                   , border  = "black"
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end if (var.plt)
      #------------------------------------------------------------------------------------#







      #------------------------------------------------------------------------------------#
      #     DBH plots by season.                                                           #
      #------------------------------------------------------------------------------------#
      if (var.plt & is.dbh){
         #---------------------------------------------------------------------------------#
         #     Loop over each season.                                                      #
         #---------------------------------------------------------------------------------#
         for (e in 1:n.season){
            cseason=paste(sprintf("%2.2i",e),season.key[e],sep="-")

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #      Plot the box plot by PFT and by season.                                 #
            #------------------------------------------------------------------------------#
            xlimit = c(0,n.pft.bp*n.simul) + 0.5
            xat    = seq(from=1 + 0.5*(n.simul-1), to=n.pft.bp*n.simul, by=n.simul)
            xgrid  = seq(from=0.5, to=(n.pft.bp+1)*n.simul - 0.5*(n.simul-1), by=n.simul)
            #------------------------------------------------------------------------------#


            #----- Load variable and set y axis. ------------------------------------------#
            this      = eft[[var.vname]]$tspftdbh[,e,,pft.mp,]
            ylimit    = pretty.xylim(u=this,fracexp=0.0,is.log=var.plog)
            #------------------------------------------------------------------------------#



            #----- Set the title. ---------------------------------------------------------#
            letitre = paste(var.desc," - ",longname,"\n"," Season: ",season.desc[e],sep="")
            lex     = paste("DBH Class [cm]")
            ley     = paste(var.desc," [",var.unit,"]",sep=" ")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Set up the box plot colour.                                              #
            #------------------------------------------------------------------------------#
            bp.colour = rep(simul$colour,times=n.dbh)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all formats.                                                  #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open file or display. -----------------------------------------------#
               fichier = paste(out.boxpftdbh[v],"/",var.vname,"-",iata,"-",cseason
                              ,"-boxplot.",outform[o],sep="")
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
               #     Split the window into several smaller windows.  Add a bottom row to   #
               # fit the legend.                                                           #
               #---------------------------------------------------------------------------#
               par.orig = par(no.readonly = TRUE)
               mar.orig = par.orig$mar
               par(oma = c(0.2,3,4,0))
               layout(mat    = rbind(1+lo.pft$mat,rep(1,times=lo.pft$ncol))
                     ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#



               #----- Plot legend. --------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
               legend ( x       = "bottom"
                      , inset   = 0.0
                      , legend  = simul$desc
                      , fill    = simul$colour
                      , border  = "black"
                      , bg      = "white"
                      , ncol    = min(3,pretty.box(n.simul)$ncol)
                      , title   = expression(bold("Simulation"))
                      , cex     = 1.0
                      )#end legend
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Loop over all seasons, and plot the bar plots.                        #
               #---------------------------------------------------------------------------#
               for (f in 1:n.pft.mp){

                  #----- Find out where is this box going, and set up axes and margins. ---#
                  left    = (f %% lo.pft$ncol) == 1
                  right   = (f %% lo.pft$ncol) == 0
                  top     = f <= lo.pft$ncol
                  bottom  = f > (lo.pft$nrow - 1) * lo.pft$ncol
                  mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
                  #------------------------------------------------------------------------#


                  #----- Set up the title for each plot. ----------------------------------#
                  lesub = paste(pft.desc[pft.mp[f]],sep="")
                  #------------------------------------------------------------------------#



                  #----- Create the temporary box plot list. ------------------------------#
                  ss.this = this[,,f,]
                  ai      = arrayInd(sequence(length(ss.this)),.dim=dim(ss.this))
                  bp.plot = split(x=ss.this,f=list(ai[,3],ai[,2]))
                  #------------------------------------------------------------------------#



                  #----- Plot window and grid. --------------------------------------------#
                  par(mar=mar.now)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
                  if (bottom) axis(side=1,at=xat,labels=dbh.key)
                  if (left  ) axis(side=2)
                  box()
                  title(main=lesub)
                  if (plotgrid) abline(v=xgrid,h=axTicks(2),col="grey40",lty="solid")

                  #----- Add the box plot, without the x axis. ----------------------------#
                  boxplot(x=bp.plot,col=bp.colour,notch=notch,add=TRUE,show.names=FALSE
                         ,yaxt="n")
                  #------------------------------------------------------------------------#
               }#end for (f in 1:n.pft.mp)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
               mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
               mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for (o in 1:nout)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         }#end for (e in 1:n.season)
         #---------------------------------------------------------------------------------#
      }#end if (var.plt && is.dbh)
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#








   #=======================================================================================#
   #=======================================================================================#
   #      Plot the bar plots.                                                              #
   #---------------------------------------------------------------------------------------#
   cat  ("   - Plotting the bar plots...","\n")
   for (v in 1:nscen.barplot){
      #----- Copy variable info. ----------------------------------------------------------#
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
      #------------------------------------------------------------------------------------#


      if (var.plt){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the bar plot by PFT, DBH class, and season.                           #
         #---------------------------------------------------------------------------------#
         cat  ("   - ",var.desc,"...","\n")
         #---------------------------------------------------------------------------------#


         #----- Load variable and set y axis. ---------------------------------------------#
         this.5d   = eft[[var.vname]]$tspftdbh[,1:n.season.mp,,pft.bp,]
         this      = apply(X=this.5d,MARGIN=c(2,3,4,5),FUN=mean,na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Stack the PFT data.                                                        #
         #---------------------------------------------------------------------------------#
         this      = apply(X=this,MARGIN=c(1,2,4),FUN=cumsum)
         this      = aperm(a=this,perm=c(4,3,1,2))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Combine simulation and size into one dimension.                            #
         #---------------------------------------------------------------------------------#
         dth  = dim(this)
         dnth = list( paste(rep(dimnames(this)[[1]],times=dth[2])
                           ,rep(dimnames(this)[[2]],each =dth[1])
                           ,sep=".")
                    , dimnames(this)[[3]]
                    , dimnames(this)[[4]]
                    )#end list
         this = array(data= this,dim=c(dth[1]*dth[2],dth[3],dth[4]),dimnames=dnth)
         #---------------------------------------------------------------------------------#


         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Seasonal means) - ",longname,sep="")
         lex     = paste("DBH class [cm]")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Find the x dimensions.                                                       #
         #---------------------------------------------------------------------------------#
         xleft  = ( rep(sequence(n.simul),times=n.dbh)
                  + rep((n.simul+1)*(sequence(n.dbh)-1),each=n.simul) ) - 0.9
         xright = xleft + 0.8
         xat    = seq(from=0.5*n.simul,by=n.simul+1,length.out=n.dbh)
         xgrid  = seq(from=0,by=n.dbh,length.out=n.simul+2)-0.5
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find out the colour filling.                                                #
         #---------------------------------------------------------------------------------#
         fill.col   = array( data     = rep(pft.colour[pft.bp],each=dim(this)[1])
                           , dim      = dim(this)[1:2]
                           , dimnames = dimnames(this)[1:2]
                           )#end array
         border.col = array( data     = rep(simul$colour,times=dim(this)[2])
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
         ybottom    = this; ybottom[,-1,] = ybottom[,-dim(ybottom)[2],]; ybottom[,1,] = 0
         ytop       = this
         #---------------------------------------------------------------------------------#




         #----- Find the limits for the plot. ---------------------------------------------#
         xlimit = range(c(xat,xgrid))
         ylimit = pretty.xylim(u=this,fracexp=0.0,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.barplot.season[v],"/",var.vname,"-",iata
                           ,"-barplot_season.",outform[o],sep="")
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
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,3,4,0))

            emat = matrix( data  = c( rep(2+t(lo.season$mat),each=2)
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
            #------------------------------------------------------------------------------#



            #----- Plot simulation legend. ------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , fill    = simul$colour
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#



            #----- Plot PFT legend. -------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = pft.desc  [pft.bp]
                   , fill    = pft.colour[pft.bp]
                   , border  = "black"
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.pft.bp)$ncol)
                   , title   = expression(bold("Plant Functional Type"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Loop over all seasons, and plot the bar plots.                           #
            #------------------------------------------------------------------------------#
            for (e in 1:n.season.mp){

               #----- Find out where is this box going, and set up axes and margins. ------#
               left    = (e %% lo.season$ncol) == 1
               right   = (e %% lo.season$ncol) == 0
               top     = e <= lo.season$ncol
               bottom  = e > (lo.season$nrow - 1) * lo.season$ncol
               mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1
               #---------------------------------------------------------------------------#


               #----- Set up the title for each plot. -------------------------------------#
               lesub = paste(season.desc[e],sep="")
               #---------------------------------------------------------------------------#



               #----- Plot window and grid. -----------------------------------------------#
               par(mar=mar.now)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
               if (bottom) axis(side=1,at=xat,labels=dbh.key)
               if (left  ) axis(side=2)
               box()
               title(main=lesub)
               if (plotgrid) abline(v=xgrid,h=axTicks(2),col="grey40",lty="solid")

               #----- Add the bar plot using rect, so we can change the line width. -------#
               rect( xleft   = xleft
                   , ybottom = ybottom [,,e]
                   , xright  = xright
                   , ytop    = ytop    [,,e]
                   , density = -1
                   , col     = fill.col
                   , border  = border.col
                   , lwd     = barplot.lwd
                   )#end rect
               #---------------------------------------------------------------------------#
            }#end for (e in 1:n.season.mp)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
            mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
         #      Plot the bar plot by PFT, DBH class, using annual means.                   #
         #---------------------------------------------------------------------------------#


         #----- Load variable and set y axis. ---------------------------------------------#
         this.4d   = eft[[var.vname]]$tspftdbh[,n.season,,pft.bp,]
         this      = apply(X=this.4d,MARGIN=c(2,3,4),FUN=mean,na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Stack the PFT data.                                                        #
         #---------------------------------------------------------------------------------#
         this      = apply(X=this,MARGIN=c(1,3),FUN=cumsum)
         this      = aperm(a=this,perm=c(3,2,1))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Combine simulation and size into one dimension.                            #
         #---------------------------------------------------------------------------------#
         dth  = dim(this)
         dnth = list( paste(rep(dimnames(this)[[1]],times=dth[2])
                           ,rep(dimnames(this)[[2]],each =dth[1])
                           ,sep=".")
                    , dimnames(this)[[3]]
                    )#end list
         this = array(data= this,dim=c(dth[1]*dth[2],dth[3]),dimnames=dnth)
         #---------------------------------------------------------------------------------#


         #----- Set the title. ------------------------------------------------------------#
         letitre = paste(var.desc," (Annual means) - ",longname,sep="")
         lex     = paste("DBH class [cm]")
         ley     = paste(var.desc," [",var.unit,"]",sep=" ")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Find the x dimensions.                                                       #
         #---------------------------------------------------------------------------------#
         xleft  = ( rep(sequence(n.simul),times=n.dbh)
                  + rep((n.simul+1)*(sequence(n.dbh)-1),each=n.simul) ) - 0.9
         xright = xleft + 0.8
         xat    = seq(from=0.5*n.simul,by=n.simul+1,length.out=n.dbh)
         xgrid  = seq(from=0,by=n.dbh,length.out=n.simul+2)-0.5
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find out the colour filling.                                                #
         #---------------------------------------------------------------------------------#
         fill.col   = array( data     = rep(pft.colour[pft.bp],each=dim(this)[1])
                           , dim      = dim(this)[1:2]
                           , dimnames = dimnames(this)[1:2]
                           )#end array
         border.col = array( data     = rep(simul$colour,times=dim(this)[2])
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
         ybottom    = this; ybottom[,-1] = ybottom[,-dim(ybottom)[2]]; ybottom[,1] = 0
         ytop       = this
         #---------------------------------------------------------------------------------#




         #----- Find the limits for the plot. ---------------------------------------------#
         xlimit = range(c(xat,xgrid))
         ylimit = pretty.xylim(u=this,fracexp=0.0,is.log=var.plog)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over all formats.                                                     #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file or display. --------------------------------------------------#
            fichier = paste(out.barplot.year[v],"/",var.vname,"-",iata
                           ,"-barplot_year.",outform[o],sep="")
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
            #     Split the window into several smaller windows.  Add a bottom row to fit  #
            # the legend.                                                                  #
            #------------------------------------------------------------------------------#
            par.orig = par(no.readonly = TRUE)
            emat = matrix(c(3,3,1,2),ncol=2,nrow=2,byrow=T)
            layout(mat=emat,heights=c(5,1))
            #------------------------------------------------------------------------------#



            #----- Plot simulation legend. ------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , fill    = simul$colour
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.simul)$ncol)
                   , title   = expression(bold("Simulation"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#



            #----- Plot PFT legend. -------------------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = pft.desc  [pft.bp]
                   , fill    = pft.colour[pft.bp]
                   , border  = "black"
                   , bg      = "white"
                   , ncol    = min(3,pretty.box(n.pft.bp)$ncol)
                   , title   = expression(bold("Plant Functional Type"))
                   , cex     = 1.0
                   )#end legend
            #------------------------------------------------------------------------------#




            #----- Plot window and grid. --------------------------------------------------#
            par(mar=mar.orig)
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit,xaxt="n",yaxt="n")
            axis(side=1,at=xat,labels=dbh.key)
            axis(side=2)
            box()
            title(main=letitre,xlab=lex,ylab=ley)
            if (plotgrid) abline(v=xgrid,h=axTicks(2),col="grey40",lty="solid")
            #------------------------------------------------------------------------------#



            #----- Add the bar plot using rect, so we can change the line width. ----------#
            rect( xleft   = xleft
                , ybottom = ybottom
                , xright  = xright
                , ytop    = ytop
                , density = -1
                , col     = fill.col
                , border  = border.col
                , lwd     = 2.5
                )#end rect
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      }#end if (var.plt)
      #------------------------------------------------------------------------------------#
   }#end for (v in 1:nscen.barplot)
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Plot the parameter space.                                                        #
   #---------------------------------------------------------------------------------------#
   cat  ("   - Plotting the parameter space...","\n")
   #---------------------------------------------------------------------------------------#
   #     Loop over Y variables.                                                            #
   #---------------------------------------------------------------------------------------#
   for (y in 1:nscen.yvar){
      y.vname = scen.xyz$yvar$vname[y]
      y.desc  = scen.xyz$yvar$desc [y]
      y.unit  = scen.xyz$yvar$unit [y]
      y.plog  = scen.xyz$yvar$plog [y]
      y.leg   = scen.xyz$yvar$leg  [y]
      y.ts    = eft[[y.vname]]$ts   [,season.mp,]
      y.tspft = eft[[y.vname]]$tspft[,n.season,pft.mp,]


      #----- Create indices to split the arrays into lists. -------------------------------#
      ai.ts    = arrayInd(sequence(length(y.ts   )),.dim=dim(y.ts   ))
      ai.tspft = arrayInd(sequence(length(y.tspft)),.dim=dim(y.tspft))
      #------------------------------------------------------------------------------------#



      #----- Find the format of the points for the plots. ---------------------------------#
      pch.ts    = array( data = simul$pch[arrayInd(1:length(y.ts   ),.dim=dim(y.ts   ))[,3]]
                       , dim  = dim(y.ts)
                       )#end array
      pch.tspft = array( data = simul$pch[arrayInd(1:length(y.tspft),.dim=dim(y.tspft))[,3]]
                       , dim  = dim(y.tspft)
                       )#end array
      #------------------------------------------------------------------------------------#




      cat("     * ",y.desc,"\n")
      #------------------------------------------------------------------------------------#
      #     Loop over X variables.                                                         #
      #------------------------------------------------------------------------------------#
      for (x in 1:nscen.xvar){
         x.vname = scen.xyz$xvar$vname[x]
         x.desc  = scen.xyz$xvar$desc [x]
         x.unit  = scen.xyz$xvar$unit [x]
         x.plog  = scen.xyz$xvar$plog [x]
         x.leg   = scen.xyz$xvar$leg  [x]
         x.ts    = eft[[x.vname]]$ts   [,season.mp,]
         x.tspft = eft[[x.vname]]$tspft[,n.season,pft.mp,]
         #---------------------------------------------------------------------------------#


         #----- Build the log info. -------------------------------------------------------#
         xy.plog=""
         if (x.plog) xy.plog=paste(xy.plog,"x",sep="")
         if (y.plog) xy.plog=paste(xy.plog,"y",sep="")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Loop over Z variables.                                                      #
         #---------------------------------------------------------------------------------#
         for (z in 1:nscen.zvar){
            z.vname   = scen.xyz$zvar$vname         [z]
            z.desc    = scen.xyz$zvar$desc          [z]
            z.unit    = scen.xyz$zvar$unit          [z]
            z.plog    = scen.xyz$zvar$plog          [z]
            z.pft     = scen.xyz$zvar$pftvar        [z]
            z.cscheme = get(scen.xyz$zvar$col.scheme[z])
            z.ts      = eft[[z.vname]]$ts   [,season.mp,]
            if (z.pft){
               z.tspft   = eft[[z.vname]]$tspft[,n.season,pft.mp,]
            }else{
               z.tspft   = array(NA,dim=dim(y.tspft),dimnames=dimnames(y.tspft))
               for (f in 1:n.pft.mp){
                  z.tspft[,f,] = eft[[z.vname]]$ts   [,n.season,]
               }#end z.tspft
            }#end 

            cat("       ~ X = ",x.desc," Z = ",z.desc,"\n")

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #     Split the arrays into lists.                                             #
            #------------------------------------------------------------------------------#
            x.list    = split(x=x.ts  ,f=ai.ts[,2])
            y.list    = split(x=y.ts  ,f=ai.ts[,2])
            z.list    = split(x=z.ts  ,f=ai.ts[,2])
            pch.list  = split(x=pch.ts,f=ai.ts[,2])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Fix the titles.                                                         #
            #------------------------------------------------------------------------------#
            letitre   = paste("Seasonal means - ",z.desc,"\n",longname,sep="")
            lex       = paste(x.desc," [",x.unit,"]",sep="")
            ley       = paste(y.desc," [",y.unit,"]",sep="")
            lacle     = paste("[",z.unit,"]",sep="")
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the parameter space by season.                                      #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open file or display. -----------------------------------------------#
               fichier = paste(out.xyz.season[y],"/y-",y.vname,"_x-",x.vname,"_z-",z.vname
                              ,"-",iata,"-season.",outform[o],sep="")
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
               #      Plot the boxes.                                                      #
               #---------------------------------------------------------------------------#
               xyz.plot( x              = x.list
                       , y              = y.list
                       , z              = z.list
                       , fixed.xlim     = TRUE
                       , fixed.ylim     = TRUE
                       , xy.log         = xy.plog
                       , pch            = pch.list
                       , cex            = 1.2
                       , nlevels        = xyz.ncolour
                       , colour.palette = z.cscheme
                       , xyz.main       = list(text=letitre)
                       , xyz.sub        = season.mp.desc
                       , xyz.xlab       = list(text=lex,adj=mtext.xadj,padj=mtext.xoff)
                       , xyz.ylab       = list(text=ley,adj=mtext.yadj,padj=mtext.yoff)
                       , xyz.more       = list(grid=list(col="grey83",lty="solid"))
                       , key.log        = z.plog
                       , key.title      = list(main=lacle,cex.main=0.8)
                       , xyz.legend     = list( x       = "bottom"
                                              , inset   = 0.0
                                              , legend  = simul$desc
                                              , pch     = simul$pch
                                              , col     = "grey16"
                                              , border  = "black"
                                              , bg      = "white"
                                              , ncol    = min(3,pretty.box(n.simul)$ncol)
                                              , title   = expression(bold("Simulation"))
                                              , cex     = 1.0
                                              )#end legend
                       )#end xyz.plot
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #     Split the arrays into lists.                                             #
            #------------------------------------------------------------------------------#
            x.list    = split(x=x.tspft  ,f=ai.tspft[,2])
            y.list    = split(x=y.tspft  ,f=ai.tspft[,2])
            z.list    = split(x=z.tspft  ,f=ai.tspft[,2])
            pch.list  = split(x=pch.tspft,f=ai.tspft[,2])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Fix the titles.                                                         #
            #------------------------------------------------------------------------------#
            letitre   = paste("Seasonal means - ",z.desc,"\n",longname,sep="")
            lex       = paste(x.desc," [",x.unit,"]",sep="")
            ley       = paste(y.desc," [",y.unit,"]",sep="")
            lacle     = paste("[",z.unit,"]",sep="")
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the parameter space by season.                                      #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open file or display. -----------------------------------------------#
               fichier = paste(out.xyz.pft[y],"/y-",y.vname,"_x-",x.vname,"_z-",z.vname
                              ,"-",iata,"-pft.",outform[o],sep="")
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
               #      Plot the boxes.                                                      #
               #---------------------------------------------------------------------------#
               xyz.plot( x              = x.list
                       , y              = y.list
                       , z              = z.list
                       , fixed.xlim     = FALSE
                       , fixed.ylim     = TRUE
                       , xy.log         = xy.plog
                       , pch            = pch.list
                       , cex            = 1.2
                       , nlevels        = xyz.ncolour
                       , colour.palette = z.cscheme
                       , xyz.main       = list(text=letitre)
                       , xyz.sub        = pft.desc[pft.mp]
                       , xyz.xlab       = list(text=lex,adj=mtext.xadj,padj=mtext.xoff)
                       , xyz.ylab       = list(text=ley,adj=mtext.yadj,padj=mtext.yoff)
                       , xyz.more       = list(grid=list(col="grey83",lty="solid"))
                       , key.log        = z.plog
                       , key.title      = list(main=lacle,cex.main=0.8)
                       , xyz.legend     = list( x       = "bottom"
                                              , inset   = 0.0
                                              , legend  = simul$desc
                                              , pch     = simul$pch
                                              , col     = "grey16"
                                              , border  = "black"
                                              , bg      = "white"
                                              , ncol    = min(3,pretty.box(n.simul)$ncol)
                                              , title   = expression(bold("Simulation"))
                                              , cex     = 1.0
                                              )#end legend
                       )#end xyz.plot
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#

}#end for (p in 1:n.sites)
#------------------------------------------------------------------------------------------#


