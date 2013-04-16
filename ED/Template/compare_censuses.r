#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#

#----- Paths. -----------------------------------------------------------------------------#
here           = getwd()                       # Current directory.
srcdir         = "/n/home00/mlongo/util/Rsc"   # Source  directory.
outdata        = file.path(here,"Rdata_census",sep="/")
#------------------------------------------------------------------------------------------#


#----- Background and output path. --------------------------------------------------------#
ibackground =     2                     # Sought background colour (actual background will
                                        #  be transparent, but foreground colours will 
                                        #  change)
                                        #  0 -- white background
                                        #  1 -- black background
                                        #  2 -- dark grey background
#----- Main path for output. --------------------------------------------------------------#
outroot     = file.path(here,paste("census_comp_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#


#----- Additional settings. ---------------------------------------------------------------#
yearend        = 2099
iphen.key      = c("iphen-01","iphen+02")
iphen.desc     = c("Evergreen","Drought deciduous")
tfall.key      = c("tfall111","tfall125","tfall140")
tfall.desc     = c("Treefall = 1.11%/yr","Treefall = 1.25%/yr","Treefall = 1.40%/yr")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = c("eps","png","pdf")   # Formats for output file.  Supported formats are:
                                        #   - "X11" - for printing on screen
                                        #   - "eps" - for postscript printing
                                        #   - "png" - for PNG printing
                                        #   - "pdf" - for PDF printing
depth          = 96                     # PNG resolution, in pixels per inch
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 16                     # Font size.
lwidth         = 2.5                    # Line width
plotgrid       = TRUE                   # Should I plot the grid in the background? 
sasfixlimits   = FALSE                  # Use a fixed scale for size and age-structure
                                        #    plots? (FALSE will set a suitable scale for
                                        #    each plot)
mtext.xoff     = -7.00                  # Offset for the x label
mtext.yoff     = -1.00                  # Offset for the y label
mtext.xadj     =  0.50                  # Offset for the x label
mtext.yadj     =  0.65                  # Offset for the y label
#------------------------------------------------------------------------------------------#


#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      =    3         # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                              # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
ed22.ci        = TRUE         # Plot confidence interval for ED?
n.boot         = 1000         # Number of realisations for bootstrap
#------------------------------------------------------------------------------------------#


#----- Simulation settings. ---------------------------------------------------------------#
place = list()
year.list = sprintf("%+2.2i",seq(from=-95,to=-5,by=10))
place[[ 1]] = list( iata   = "gyf" 
                  , config = list( yeara    = paste("yra",year.list,sep="")
                                 , iphen    = iphen.key
                                 , stext    = c("stext06","stext08")
                                 , treefall = tfall.key
                                 )#end list
                  )#end list
place[[ 2]] = list( iata   = "s67" 
                  , config = list( yeara    = paste("yra",year.list,sep="")
                                 , iphen    = iphen.key
                                 , stext    = c("stext16","stext11")
                                 , treefall = tfall.key
                                 )#end list
                  )#end list
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#



#----- Loading some packages and scripts. -------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Types of variables to use to determine mortality, growth, and recruitment.            #
#------------------------------------------------------------------------------------------#
recr.vars     = c("n","agb","ba")
recr.labels   = c("Individuals","Above Ground Biomass","Basal area")
recr.units    = c(untab$pcpopoyr,untab$pcagboyr,untab$pcbaoyr)
mort.vars     = c("n","agb","ba")
mort.labels   = c("Individuals","Above Ground Biomass","Basal area")
mort.units    = c(untab$pcpopoyr,untab$pcagboyr,untab$pcbaoyr)
growth.vars   = c("dbh","agb","ba")
growth.labels = c("DBH","Above Ground Biomass","Basal Area")
growth.units  = c(untab$pcdbhoyr,untab$pcagboyr,untab$pcbaoyr)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Comparisons.                                                                         #
#------------------------------------------------------------------------------------------#
#---- 1. Plot time series of expected values and confidence intervals. --------------------#
pratets      = list()
pratets[[1]] = list( ed2.rate   = "recr"
                   , sta.rate   = "recr"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Recruitment rate"
                   , unit.rate  = recr.units
                   , col.ed2    = c(chartreuse.fg,chartreuse.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = recr.vars
                   , desc.indiv = recr.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[2]] = list( ed2.rate   = "mort"
                   , sta.rate   = "mort"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Mortality rate"
                   , unit.rate  = mort.units
                   , col.ed2    = c(purple.fg,purple.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = mort.vars
                   , desc.indiv = mort.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[3]] = list( ed2.rate   = "ddmort"
                   , sta.rate   = "mort"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Density-dependent mort. rate"
                   , unit.rate  = mort.units
                   , col.ed2    = c(indigo.fg,indigo.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = mort.vars
                   , desc.indiv = mort.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[4]] = list( ed2.rate   = "dimort"
                   , sta.rate   = "mort"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Density-independent mort. rate"
                   , unit.rate  = mort.units
                   , col.ed2    = c(blue.fg,blue.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = mort.vars
                   , desc.indiv = mort.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[5]] = list( ed2.rate   = "growth"
                   , sta.rate   = "growth"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Growth rate"
                   , unit.rate  = growth.units
                   , col.ed2    = c(yellow.fg,yellow.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = growth.vars
                   , desc.indiv = growth.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list



#---- 2. Plot expected values and confidence intervals for all size classes and censuses. -#
pratesize      = list()
pratesize[[1]] = list( ed2.rate   = "recr"
                     , sta.rate   = "recr"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Recruitment rate"
                     , unit.rate  = recr.units
                     , col.ed2    = c(chartreuse.fg,chartreuse.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = recr.vars
                     , desc.indiv = recr.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[2]] = list( ed2.rate   = "mort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Mortality rate"
                     , unit.rate  = mort.units
                     , col.ed2    = c(purple.fg,purple.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = mort.vars
                     , desc.indiv = mort.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[3]] = list( ed2.rate   = "ddmort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Density-dependent mort. rate"
                     , unit.rate  = mort.units
                     , col.ed2    = c(indigo.fg,indigo.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = mort.vars
                     , desc.indiv = mort.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[4]] = list( ed2.rate   = "dimort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Density-independent mort. rate"
                     , unit.rate  = mort.units
                     , col.ed2    = c(blue.fg,blue.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = mort.vars
                     , desc.indiv = mort.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[5]] = list( ed2.rate   = "growth"
                     , sta.rate   = "growth"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Growth rate"
                     , unit.rate  = growth.units
                     , col.ed2    = c(yellow.fg,yellow.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = growth.vars
                     , desc.indiv = growth.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
#---- 3. Plot expected values and confidence intervals for themes. ------------------------#
pratetheme      = list()
pratetheme[[1]] = list( ed2.rate   = c("ddmort","dimort")
                      , sta.rate   = "mort"
                      , sizetoo    = TRUE
                      , pfttoo     = TRUE
                      , desc.rate  = c("ED-2.2 Negative C Balance"
                                      ,"ED-2.2 Other")
                      , unit.rate  = mort.units
                      , col.ed2    = rbind( c(indigo.fg,indigo.bg)
                                          , c(chartreuse.fg ,chartreuse.bg)
                                          )#end rbind
                      , col.sta    = c(grey.fg,grey.bg)
                      , angle      = c(-45,45)
                      , density    = c(40,40)
                      , indiv      = mort.vars
                      , desc.indiv = mort.labels
                      , theme      = "mortality"
                      , theme.desc = "Mortality Rates"
                      , plog       = ""
                      )#end list
#------------------------------------------------------------------------------------------#



#----- Find how many treefall and phenology runs were set. --------------------------------#
n.tfall = length(tfall.key)
n.iphen = length(iphen.key)
#------------------------------------------------------------------------------------------#



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#


#----- Set how many variables we will compare. --------------------------------------------#
npratesize  = length(pratesize )
npratetheme = length(pratetheme)
npratets    = length(pratets   )
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#


#----- Load census data. ------------------------------------------------------------------#
census.file = paste(srcdir,"LBA_MIP.census_summ.RData",sep="/")
load(file=census.file)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Types of variables to use to determine mortality, growth, and recruitment.            #
#------------------------------------------------------------------------------------------#
nrecr.vars   = length(recr.vars  ) 
nmort.vars   = length(mort.vars  ) 
ngrowth.vars = length(growth.vars) 
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size            = plotsize(proje=FALSE,paper=paper)
wide.size       = size
wide.size$width = 1.33 * size$width
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#


#----- List of places. --------------------------------------------------------------------#
nplaces = length(place)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (p in sequence(nplaces)){
   #----- Build the list of site runs. ----------------------------------------------------#
   iata      = place[[p]]$iata
   idx       = match(iata,poilist$iata)
   longname  = poilist$longname[idx]
   config    = expand.grid(place[[p]]$config,stringsAsFactors=FALSE)
   simul.all = paste("t",iata,"_",apply(X=config,MARGIN=1,FUN=paste,collapse="_")
                    ,sep="")



   cat (" + Site: ",longname,"\n")

   #---------------------------------------------------------------------------------------#
   #     Treefall loop.                                                                    #
   #---------------------------------------------------------------------------------------#
   for (tf in sequence(n.tfall)){
      #----- Set up the path for this treefall. -------------------------------------------#
      outfall = file.path(outroot,tfall.key[tf])
      if (! file.exists(outfall)) dir.create(outfall)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Phenology loop.                                                                #
      #------------------------------------------------------------------------------------#
      for (ph in sequence(n.iphen)){

         cat ("  - Group: ",tfall.desc[tf]," - ",iphen.desc[ph],"\n")



         #----- Use only this phenology and treefall combination. -------------------------#
         aux.iphen.key=sub(pattern="\\+",replacement=".",x=iphen.key)
         aux.tfall.key=sub(pattern="\\+",replacement=".",x=tfall.key)
         aux.simul.all=sub(pattern="\\+",replacement=".",x=simul.all)
         keep    = ( regexpr(pattern=aux.iphen.key[ph],text=aux.simul.all) > 0
                   & regexpr(pattern=aux.tfall.key[tf],text=aux.simul.all) > 0 )
         simul   = simul.all[keep]
         n.simul = length(simul)
         #---------------------------------------------------------------------------------#


         #----- Set up the location. ------------------------------------------------------#
         lieu    = paste(longname,iphen.desc[ph],tfall.desc[tf],sep=" - ")
         #---------------------------------------------------------------------------------#




         #----- Find the census observations for this particular site. --------------------#
         if (iata == "mao" | iata == "bdf"){
            census.name = "census.m34"
         }else if(iata == "stm" | iata == "s66"){
            census.name = "census.s67"
         }else if(iata == "rao"){
            census.name = "census.pdg"
         }else if(iata == "jpr"){
            census.name = "census.fns"
         }else if(iata == "btr"){
            census.name = "census.s77"
         }else{
            census.name = paste("census.",iata,sep="")
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     We only run this part of the code if there are observations to compare with #
         # the model.                                                                      #
         #---------------------------------------------------------------------------------#
         if (census.name %in% ls()){

            #------------------------------------------------------------------------------#
            #     Load the census data, from the monthly means.                            #
            #------------------------------------------------------------------------------#
            sta         = get(census.name)
            n.census    = length(sta$when)
            n.dbh       = length(sta$dbh.breaks)-1
            
            x.edge      = c(10,      sta$dbh.breaks[-c(1,n.dbh+1)]
                              , 2. * sta$dbh.breaks[n.dbh] - sta$dbh.breaks[n.dbh-1] )
            x.dbh       = 0.5 * ( x.edge[-1] + x.edge[-(n.dbh+1)] )
            xlimit      = pretty.xylim(u=x.edge,fracexp=0.0,is.log=FALSE)
            dbh.names   = dimnames(sta$mort.size$n$expected)[[2]]
            year4       = numyears(sta$when)
            biocyca     = year4[2]
            biocycz     = year4[length(year4)]
            dyear       = c(NA,diff(year4))
            census.desc = paste(year4-c(NA,diff(year4)),year4,sep="-")


            #------------------------------------------------------------------------------#
            #      Loop over all months to grab all the census data.                       #
            #------------------------------------------------------------------------------#
            census.idx   = NULL
            for (y in 2:n.census){
               #----- Find the first and last time to be averaged for this census. --------#
               ts.montha  = ( nummonths(sta$when[y-1]) %% 12 )
               ts.yeara   = numyears (sta$when[y-1])
               ts.monthz  = ( ( (nummonths(sta$when[y]) - 1) %% 12 )
                            + 12 * as.integer(nummonths(sta$when[y]) == 1) )
               ts.yearz   = numyears (sta$when[y]) - as.integer(ts.monthz == 12)
               n.inter    = (ts.yearz-ts.yeara-1)*12 + ts.monthz + (12 - ts.montha + 1)
               #---------------------------------------------------------------------------#
               census.idx = c(census.idx,rep(y,times=n.inter))
            }#end for
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop over all simulations to build the data frames.                      #
            #------------------------------------------------------------------------------#
            for (s in sequence(n.simul)){
               cat ("      * Get data from simulation ",simul[s],"...","\n")


               #---------------------------------------------------------------------------#
               #     Retrieve default information about this place and set up some vari-   #
               # ables.                                                                    #
               #---------------------------------------------------------------------------#
               yeara   = as.numeric(substring(config$yeara[s],4,7))
               thispoi = locations(where=simul[s],here=here,yearbeg=yeara,yearend=yearend
                                  ,monthbeg=1)
               inpref  = thispoi$pathin
               outmain = file.path(here,simul[s])
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Load the data for this simulation.                                    #
               #---------------------------------------------------------------------------#
               this.rdata = file.path(outmain,"rdata_census"
                                     ,paste("census_",simul[s],".RData",sep=""))
               load(this.rdata)
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     If this is the first time, create the temporary vectors.              #
               #---------------------------------------------------------------------------#
               if (s == 1){
                  n.months = ed2$tseries$n.months
                  n.cycles = ed2$tseries$n.cycles

                  #------------------------------------------------------------------------#
                  #     Make the vector with the environmental, plot-level, and size-level #
                  # dimensions.                                                            #
                  #------------------------------------------------------------------------#
                  dim.size     = c(npft+1,n.dbh+1,n.months,n.cycles,n.simul)
                  #----- Recruitment. -----------------------------------------------------#
                  ts.recr.size     = list()
                  for (v in 1:nrecr.vars){
                     ts.recr.size  [[recr.vars[v]]] = array( NA, dim = dim.size)
                  }#end for
                  #----- Mortality. -------------------------------------------------------#
                  ts.mort.size     = list()
                  ts.ddmort.size   = list()
                  ts.dimort.size   = list()
                  for (v in 1:nmort.vars){
                     ts.mort.size  [[mort.vars[v]]] = array( NA, dim = dim.size)
                     ts.ddmort.size[[mort.vars[v]]] = array( NA, dim = dim.size)
                     ts.dimort.size[[mort.vars[v]]] = array( NA, dim = dim.size)
                  }#end for
                  #----- Growth. ----------------------------------------------------------#
                  ts.growth.size   = list()
                  for (v in 1:ngrowth.vars){
                     ts.growth.size[[growth.vars[v]]] = array( NA, dim = dim.size)
                  }#end for
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#



               #----- Copy the information from this simulation. --------------------------#
               for (v in 1:nrecr.vars){
                  var.now = recr.vars[v]
                  ts.recr.size  [[var.now]][,,,,s] = ed2$tseries$recr  [[var.now]]
               }#end for
               for (v in 1:nmort.vars){
                  var.now = mort.vars[v]
                  ts.mort.size  [[var.now]][,,,,s] = ed2$tseries$mort  [[var.now]]
                  ts.ddmort.size[[var.now]][,,,,s] = ed2$tseries$ddmort[[var.now]]
                  ts.dimort.size[[var.now]][,,,,s] = ed2$tseries$dimort[[var.now]]
               }#end for
               for (v in 1:ngrowth.vars){
                  var.now = growth.vars[v]
                  ts.growth.size[[var.now]][,,,,s] = ed2$tseries$growth[[var.now]]
               }#end for
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     Delete the simulation because ed2 becomes the combination of all      #
               # simulations.                                                              #
               #---------------------------------------------------------------------------#
               rm(ed2)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(n.simul)
            #------------------------------------------------------------------------------#






            #==============================================================================#
            #==============================================================================#
            #------------------------------------------------------------------------------#
            #     Remove the 4th and 5th. dimensions for the means.                        #
            #------------------------------------------------------------------------------#
            cat("   - Averaging the census intervals...","\n")


            #----- This function will be used by the bootstrap. ---------------------------#
            mean.fun = function(x,idx) mean(x[idx],na.rm=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot-level recruitment.                                                  #
            #------------------------------------------------------------------------------#
            ms.recr.plot   = list()
            ms.mort.plot   = list()
            ms.ddmort.plot = list()
            ms.dimort.plot = list()
            ms.growth.plot = list()
            ms.recr.size   = list()
            ms.mort.size   = list()
            ms.ddmort.size = list()
            ms.dimort.size = list()
            ms.growth.size = list()
            #------------------------------------------------------------------------------#



            #----- Plot-level recruitment. ------------------------------------------------#
            cat("     * Recruitment...","\n")
            for (v in 1:nrecr.vars){
               v.now  = recr.vars[v]
               ms.now = list()

               recr.rates = c("recr")
               n.rates    = length(recr.rates)

               #---------------------------------------------------------------------------#
               #     Loop over all rates.                                                  #
               #---------------------------------------------------------------------------#
               for (r in 1:n.rates){
                  ts.this.size = get(paste("ts",recr.rates[r],"size",sep="."))
                  ms.this.plot = get(paste("ms",recr.rates[r],"plot",sep="."))
                  ms.this.size = get(paste("ms",recr.rates[r],"size",sep="."))


                  #------------------------------------------------------------------------#
                  #     Find the mean for all PFTs.                                        #
                  #------------------------------------------------------------------------#
                  ms.mean.plot = array(NA,dim=c(npft+1,n.census))
                  ms.q025.plot = array(NA,dim=c(npft+1,n.census))
                  ms.q975.plot = array(NA,dim=c(npft+1,n.census))
                  ms.mean.size = array(NA,dim=c(npft+1,n.dbh,n.census))
                  ms.q025.size = array(NA,dim=c(npft+1,n.dbh,n.census))
                  ms.q975.size = array(NA,dim=c(npft+1,n.dbh,n.census))

                  for (p in 1:(npft+1)){
                     for (i in 2:n.census){
                        i.sel = census.idx == i
                        ts.plot.now       = c(ts.this.size[[v.now]][p,ndbh+1,i.sel,,])
                        ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                        if (any(is.finite(ts.plot.now))){
                           boot.now = boot   (data=ts.plot.now,statistic=mean.fun,R=n.boot)
                           ci.now   = try(boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                                         ,silent=TRUE)
                           if ("try-error" %in% is(ci.now)){
                              warning("Failed using bootstrap...")
                           }else if (length(ci.now$percent) == 5){
                              ms.q025.plot[p,i] = ci.now$percent[4]
                              ms.q975.plot[p,i] = ci.now$percent[5]
                           }else{
                              warning("Failed using bootstrap...")
                           }#end if
                        }#end if
                        #------------------------------------------------------------------#

                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           ts.size.now         = c(ts.this.size[[v.now]][p,d,i.sel,,])
                           ms.mean.size[p,d,i] = mean(ts.size.now,na.rm=TRUE) 
                           if (any(is.finite(ts.size.now))){
                              boot.now = boot(data=ts.size.now,statistic=mean.fun,R=n.boot)
                              ci.now   = try(boot.ci(boot.out=boot.now,conf=0.95
                                                    ,type="perc")
                                            ,silent=TRUE)
                              if ("try-error" %in% is(ci.now)){
                                 warning("Failed using bootstrap...")
                              }else if (length(ci.now$percent) == 5){
                                 ms.q025.size[p,d,i] = ci.now$percent[4]
                                 ms.q975.size[p,d,i] = ci.now$percent[5]
                              }else{
                                 warning("Failed using bootstrap...")
                              }#end if
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for
                        #------------------------------------------------------------------#
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#



                  #----- Convert rates fraction rates. ------------------------------------#
                  ms.mean.plot = exp(ms.mean.plot) - 1.
                  ms.q025.plot = exp(ms.q025.plot) - 1.
                  ms.q975.plot = exp(ms.q975.plot) - 1.
                  ms.mean.size = exp(ms.mean.size) - 1.
                  ms.q025.size = exp(ms.q025.size) - 1.
                  ms.q975.size = exp(ms.q975.size) - 1.
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Save plot and size.                                                #
                  #------------------------------------------------------------------------#
                  ms.this.plot[[v.now]] = list( expected = ms.mean.plot
                                              , q025     = ms.q025.plot
                                              , q975     = ms.q975.plot
                                              )#end list
                  ms.this.size[[v.now]] = list( expected = ms.mean.size
                                              , q025     = ms.q025.size
                                              , q975     = ms.q975.size
                                              )#end list
                  #------------------------------------------------------------------------#
                  dummy = assign(paste("ms",recr.rates[r],"plot",sep="."), ms.this.plot)
                  dummy = assign(paste("ms",recr.rates[r],"size",sep="."), ms.this.size)
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#


            #----- Mortality rates. -------------------------------------------------------#
            cat("     * Mortality...","\n")
            for (v in 1:nmort.vars){
               v.now  = mort.vars[v]
               ms.now = list()

               mort.rates = c("mort","ddmort","dimort")
               n.rates    = length(mort.rates)

               #---------------------------------------------------------------------------#
               #     Loop over all rates.                                                  #
               #---------------------------------------------------------------------------#
               for (r in 1:n.rates){
                  ts.this.size = get(paste("ts",mort.rates[r],"size",sep="."))
                  ms.this.plot = get(paste("ms",mort.rates[r],"plot",sep="."))
                  ms.this.size = get(paste("ms",mort.rates[r],"size",sep="."))


                  #------------------------------------------------------------------------#
                  #     Find the mean for all PFTs.                                        #
                  #------------------------------------------------------------------------#
                  ms.mean.plot = array(NA,dim=c(npft+1,n.census))
                  ms.q025.plot = array(NA,dim=c(npft+1,n.census))
                  ms.q975.plot = array(NA,dim=c(npft+1,n.census))
                  ms.mean.size = array(NA,dim=c(npft+1,n.dbh,n.census))
                  ms.q025.size = array(NA,dim=c(npft+1,n.dbh,n.census))
                  ms.q975.size = array(NA,dim=c(npft+1,n.dbh,n.census))

                  for (p in 1:(npft+1)){
                     for (i in 2:n.census){
                        i.sel = census.idx == i
                        ts.plot.now       = c(ts.this.size[[v.now]][p,ndbh+1,i.sel,,])
                        ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                        if (any(is.finite(ts.plot.now))){
                           boot.now = boot   (data=ts.plot.now,statistic=mean.fun,R=n.boot)
                           ci.now   = try(boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                                         ,silent=TRUE)
                           if ("try-error" %in% is(ci.now)){
                              warning("Failed using bootstrap...")
                           }else if (length(ci.now$percent) == 5){
                              ms.q025.plot[p,i] = ci.now$percent[4]
                              ms.q975.plot[p,i] = ci.now$percent[5]
                           }else{
                              warning("Failed using bootstrap...")
                           }#end if
                        }#end if
                        #------------------------------------------------------------------#

                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           ts.size.now         = c(ts.this.size[[v.now]][p,d,i.sel,,])
                           ms.mean.size[p,d,i] = mean(ts.size.now,na.rm=TRUE) 
                           if (any(is.finite(ts.size.now))){
                              boot.now = boot(data=ts.size.now,statistic=mean.fun,R=n.boot)
                              ci.now   = try(boot.ci(boot.out=boot.now,conf=0.95
                                                    ,type="perc")
                                            ,silent=TRUE)
                              if ("try-error" %in% is(ci.now)){
                                 warning("Failed using bootstrap...")
                              }else if (length(ci.now$percent) == 5){
                                 ms.q025.size[p,d,i] = ci.now$percent[4]
                                 ms.q975.size[p,d,i] = ci.now$percent[5]
                              }else{
                                 warning("Failed using bootstrap...")
                              }#end if
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for
                        #------------------------------------------------------------------#
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#



                  #----- Convert rates fraction rates. ------------------------------------#
                  ms.mean.plot = 1. - exp( - ms.mean.plot)
                  ms.q025.plot = 1. - exp( - ms.q025.plot)
                  ms.q975.plot = 1. - exp( - ms.q975.plot)
                  ms.mean.size = 1. - exp( - ms.mean.size)
                  ms.q025.size = 1. - exp( - ms.q025.size)
                  ms.q975.size = 1. - exp( - ms.q975.size)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Save plot and size.                                                #
                  #------------------------------------------------------------------------#
                  ms.this.plot[[v.now]] = list( expected = ms.mean.plot
                                              , q025     = ms.q025.plot
                                              , q975     = ms.q975.plot
                                              )#end list
                  ms.this.size[[v.now]] = list( expected = ms.mean.size
                                              , q025     = ms.q025.size
                                              , q975     = ms.q975.size
                                              )#end list
                  #------------------------------------------------------------------------#
                  dummy = assign(paste("ms",mort.rates[r],"plot",sep="."), ms.this.plot)
                  dummy = assign(paste("ms",mort.rates[r],"size",sep="."), ms.this.size)
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#

            
            #----- Growth rates. ----------------------------------------------------------#
            cat("     * Growth...","\n")
            for (v in 1:ngrowth.vars){
               v.now  = growth.vars[v]
               ms.now = list()

               growth.rates = c("growth")
               n.rates    = length(growth.rates)

               #---------------------------------------------------------------------------#
               #     Loop over all rates.                                                  #
               #---------------------------------------------------------------------------#
               for (r in 1:n.rates){
                  ts.this.size = get(paste("ts",growth.rates[r],"size",sep="."))
                  ms.this.plot = get(paste("ms",growth.rates[r],"plot",sep="."))
                  ms.this.size = get(paste("ms",growth.rates[r],"size",sep="."))

                  #------------------------------------------------------------------------#
                  #     Find the mean for all PFTs.                                        #
                  #------------------------------------------------------------------------#
                  ms.mean.plot = array(NA,dim=c(npft+1,n.census))
                  ms.q025.plot = array(NA,dim=c(npft+1,n.census))
                  ms.q975.plot = array(NA,dim=c(npft+1,n.census))
                  ms.mean.size = array(NA,dim=c(npft+1,n.dbh,n.census))
                  ms.q025.size = array(NA,dim=c(npft+1,n.dbh,n.census))
                  ms.q975.size = array(NA,dim=c(npft+1,n.dbh,n.census))

                  for (p in 1:(npft+1)){
                     for (i in 2:n.census){
                        i.sel = census.idx == i
                        ts.plot.now       = c(ts.this.size[[v.now]][p,ndbh+1,i.sel,,])
                        ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                        if (any(is.finite(ts.plot.now))){
                           boot.now = boot   (data=ts.plot.now,statistic=mean.fun,R=n.boot)
                           ci.now   = try(boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                                         ,silent=TRUE)
                           if ("try-error" %in% is(ci.now)){
                              warning("Failed using bootstrap...")
                           }else if (length(ci.now$percent) == 5){
                              ms.q025.plot[p,i] = ci.now$percent[4]
                              ms.q975.plot[p,i] = ci.now$percent[5]
                           }else{
                              warning("Failed using bootstrap...")
                           }#end if
                        }#end if
                        #------------------------------------------------------------------#

                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           ts.size.now         = c(ts.this.size[[v.now]][p,d,i.sel,,])
                           ms.mean.size[p,d,i] = mean(ts.size.now,na.rm=TRUE) 
                           if (any(is.finite(ts.size.now))){
                              boot.now = boot(data=ts.size.now,statistic=mean.fun,R=n.boot)
                              ci.now   = try(boot.ci(boot.out=boot.now,conf=0.95
                                                    ,type="perc")
                                            ,silent=TRUE)
                              if ("try-error" %in% is(ci.now)){
                                 warning("Failed using bootstrap...")
                              }else if (length(ci.now$percent) == 5){
                                 ms.q025.size[p,d,i] = ci.now$percent[4]
                                 ms.q975.size[p,d,i] = ci.now$percent[5]
                              }else{
                                 warning("Failed using bootstrap...")
                              }#end if
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for
                        #------------------------------------------------------------------#
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Save plot and size.                                                #
                  #------------------------------------------------------------------------#
                  ms.this.plot[[v.now]] = list( expected = ms.mean.plot
                                              , q025     = ms.q025.plot
                                              , q975     = ms.q975.plot
                                              )#end list
                  ms.this.size[[v.now]] = list( expected = ms.mean.size
                                              , q025     = ms.q025.size
                                              , q975     = ms.q975.size
                                              )#end list
                  #------------------------------------------------------------------------#
                  dummy = assign(paste("ms",growth.rates[r],"plot",sep="."), ms.this.plot)
                  dummy = assign(paste("ms",growth.rates[r],"size",sep="."), ms.this.size)
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Find the average rates for the census period using log-normal.           #
            #------------------------------------------------------------------------------#
            mypfts = sort(match(unique(names(sta$classes)),pft$name))
            npfts  = length(mypfts)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Retrieve the factor, classes and wood density used for the obser-       #
            # vations.  We can't switch the factors between observation and statistics     #
            # because we use observations to drive the statistics.                         #
            #------------------------------------------------------------------------------#
            cat("   - Finding the average rates...","\n")
            ed2           = list()
            ed2$when      = sta$when
            ed2$taxon     = sta$taxon
            ed2$classes   = sta$classes
            nfac          = length(ed2$classes)
            ed2$wood.dens = sta$wood.dens
            ed2$dtime     = sta$dtime
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Recruitment rates.  Only global values can be found.                    #
            #------------------------------------------------------------------------------#
            for (r in 1:npratets){
               #---------------------------------------------------------------------------#
               #     Load the rate information.                                            #
               #---------------------------------------------------------------------------#
               this.rate = pratets[[ r]]
               ed2.rate  = this.rate$ed2.rate
               sta.rate  = this.rate$sta.rate
               sizetoo   = this.rate$sizetoo
               desc.rate = this.rate$desc.rate
               indiv     = this.rate$indiv
               cat(" - Compounding the ",desc.rate," tables...","\n")
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Load plot-level and size-level data.                                 #
               #---------------------------------------------------------------------------#
               ed2.plot   = paste(ed2.rate,"plot",sep=".")
               sta.plot   = paste(sta.rate,"plot",sep=".")
               ed2.size   = paste(ed2.rate,"size",sep=".")
               sta.size   = paste(sta.rate,"size",sep=".")
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #       Find how many individuals to retrieve.                              #
               #---------------------------------------------------------------------------#
               nindiv = length(indiv)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Loop over the different types of individuals.                         #
               #---------------------------------------------------------------------------#
               ed2[[ed2.plot]] = list()
               for (v in 1:nindiv){
                  #----- Set up the individuals. ------------------------------------------#
                  vn   = indiv[v]
                  yy   = 2:n.census
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot-level with all cycles.                                        #
                  #------------------------------------------------------------------------#
                  ms.plot   = get(paste("ms",ed2.rate,"plot",sep="."))[[vn]]
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Create the plot-level structure.                                   #
                  #------------------------------------------------------------------------#
                  ed2[[ed2.plot]][[vn]]               = list()
                  #----- "Borrow" the structure from the sta counterpart. -----------------#
                  ed2[[ed2.plot]][[vn]]$global   = NA * sta[[sta.plot]][[vn]]$global
                  ed2[[ed2.plot]][[vn]]$expected = NA * sta[[sta.plot]][[vn]]$expected
                  ed2[[ed2.plot]][[vn]]$q025     = NA * sta[[sta.plot]][[vn]]$q025
                  ed2[[ed2.plot]][[vn]]$q975     = NA * sta[[sta.plot]][[vn]]$q975
                  #----- Save the global variables. ---------------------------------------#
                  ed2[[ed2.plot]][[vn]]$global  [1,yy] = ms.plot$expected[npft+1,yy]
                  ed2[[ed2.plot]][[vn]]$global  [2,yy] = ms.plot$q025    [npft+1,yy]
                  ed2[[ed2.plot]][[vn]]$global  [3,yy] = ms.plot$q975    [npft+1,yy]
                  #----- Save the PFT statistics. -----------------------------------------#
                  ed2[[ed2.plot]][[vn]]$expected[ ,yy] = ms.plot$expected[mypfts,yy]
                  ed2[[ed2.plot]][[vn]]$q025    [ ,yy] = ms.plot$q025    [mypfts,yy]
                  ed2[[ed2.plot]][[vn]]$q975    [ ,yy] = ms.plot$q975    [mypfts,yy]
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #    Now the size-dependent variables, if this is a size-dependent rate. #
                  #------------------------------------------------------------------------#
                  if (sizetoo){

                     #---------------------------------------------------------------------#
                     #     Plot-level with all cycles.                                     #
                     #---------------------------------------------------------------------#
                     ms.size   = get(paste("ms",ed2.rate,"size",sep="."))[[vn]]
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Create the plot-level structure.                                #
                     #---------------------------------------------------------------------#
                     ed2[[ed2.size]][[vn]]        = list()
                     #----- "Borrow" the structure from the sta counterpart. --------------#
                     ed2[[ed2.size]][[vn]]$global   = NA * sta[[sta.size]][[vn]]$global
                     ed2[[ed2.size]][[vn]]$expected = NA * sta[[sta.size]][[vn]]$expected
                     ed2[[ed2.size]][[vn]]$q025     = NA * sta[[sta.size]][[vn]]$q025
                     ed2[[ed2.size]][[vn]]$q975     = NA * sta[[sta.size]][[vn]]$q975
                     #----- Save the global variables. ------------------------------------#
                     ed2[[ed2.size]][[vn]]$global[1,,yy] = ms.size$expected[npft+1,,yy]
                     ed2[[ed2.size]][[vn]]$global[2,,yy] = ms.size$q025    [npft+1,,yy]
                     ed2[[ed2.size]][[vn]]$global[3,,yy] = ms.size$q975    [npft+1,,yy]
                     #----- Save the PFT statistics. --------------------------------------#
                     ed2[[ed2.size]][[vn]]$expected[ ,,yy] = ms.size$expected[mypfts,,yy]
                     ed2[[ed2.size]][[vn]]$q025    [ ,,yy] = ms.size$q025    [mypfts,,yy]
                     ed2[[ed2.size]][[vn]]$q975    [ ,,yy] = ms.size$q975    [mypfts,,yy]
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the directories.                                                    #
            #------------------------------------------------------------------------------#
            outplot   = paste(outfall,"census_plot"  ,sep="/")
            outsize   = paste(outfall,"census_size"  ,sep="/")
            if (! file.exists(outplot)) dir.create(outplot)
            if (! file.exists(outsize)) dir.create(outsize)
            #------------------------------------------------------------------------------#






            #==============================================================================#
            #==============================================================================#
            #  6.  Plot rates as function of size.                                         #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratesize){
               this.plot  = pratesize[[n]]
               desc.rate  = this.plot$desc.rate
               unit.rate  = this.plot$unit.rate
               col.ed2    = this.plot$col.ed2
               col.sta    = this.plot$col.sta
               indiv      = this.plot$indiv
               desc.indiv = this.plot$desc.indiv
               legpos     = this.plot$legpos
               plog       = this.plot$plog
               ylog       = length(grep("y",plog)) > 0

               nindiv     = length(indiv)



               #----- Create a directory for this type of plot. ---------------------------#
               outrate   = paste(outsize,ed2.rate,sep="/")
               if (! file.exists(outrate)) dir.create(outrate)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    Loop over all possible types of population count.                      #
               #---------------------------------------------------------------------------#
               for (i in 1:nindiv){
                  cat("  - Size level: ",desc.indiv[i],"...","\n")


                  #----- Build the rate name. ---------------------------------------------#
                  cat(" + Plotting size-dependent ",desc.rate,"...","\n")
                  ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                  sta.rate = paste(this.plot$sta.rate,"size",sep=".")
                  #------------------------------------------------------------------------#



                  #----- Create path for this individual. ---------------------------------#
                  outindiv   = paste(outrate,indiv[i],sep="/")
                  if (! file.exists(outindiv)) dir.create(outindiv)
                  #------------------------------------------------------------------------#


                  #----- Load the modelled rates. -----------------------------------------#
                  sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                  sta.expected = 100. * sta.mod[1,,]
                  sta.q025     = 100. * sta.mod[2,,]
                  sta.q975     = 100. * sta.mod[3,,]
                  ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
                  ed2.expected = 100. * ed2.mod[1,,]
                  ed2.q025     = 100. * ed2.mod[2,,]
                  ed2.q975     = 100. * ed2.mod[3,,]
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Define a nice configuration for the multiple panels.              #
                  #------------------------------------------------------------------------#
                  lo.box = pretty.box(n=n.census-1,horizontal=TRUE)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #      Loop over all formats.                                            #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Open the file or the plot window. -----------------------------#
                     fichier = file.path(outindiv
                                        ,paste(iata,"-yrsize-",ed2.rate,"-",indiv[i]
                                              ,"-",iphen.key[ph],".",outform[o] ,sep=""))
                     if(outform[o] == "x11"){
                        X11(width=wide.size$width,height=wide.size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=wide.size$width*depth
                           ,height=wide.size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=wide.size$width
                                  ,height=wide.size$height,pointsize=ptsz
                                  ,paper=wide.size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE
                           ,width=wide.size$width,height=wide.size$height,pointsize=ptsz
                           ,paper=wide.size$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Split the window into several smaller windows.                  #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     par.orig = par(no.readonly = TRUE)
                     par(oma = c(0.2,3,4,0))
                     layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                           ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                           )#end layout
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Find the plot limit for the y scale.                            #
                     #---------------------------------------------------------------------#
                     if (ed22.ci){
                        yuse   = c(sta.q025,sta.q975,ed2.q025,ed2.q975)
                        ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)



                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0.01
                               , legend  = c("Census","ED-2.2")
                               , fill    = c(col.sta[2],col.ed2[2])
                               , border  = c(col.sta[2],col.ed2[2])
                               , col     = c(col.sta[1],col.ed2[1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(-45,45)
                               , density = c( 40,40)
                               , bg      = background
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 1.0
                               )#end legend
                        #------------------------------------------------------------------#
                     }else{
                        yuse   = c(sta.q025,sta.q975,ed2.expected)
                        ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)



                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0.01
                               , legend  = c("Census","ED-2.2")
                               , fill    = c(col.sta[2],         0)
                               , border  = c(col.sta[2],         0)
                               , col     = c(col.sta[1],col.ed2[1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(-45,45)
                               , density = c( 40,40)
                               , bg      = background
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 1.0
                               )#end legend
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Loop over all years.                                           #
                     #---------------------------------------------------------------------#
                     for (y in 2:n.census){
                        k       = y - 1
                        left    = (k %% lo.box$ncol) == 1
                        right   = (k %% lo.box$ncol) == 0
                        top     = k <= lo.box$ncol
                        bottom  = k > (lo.box$nrow - 1) * lo.box$ncol
                        mar.now = c(1.1 + 4 * bottom
                                   ,1.1 + 3 * left
                                   ,1.1 + 4 * top
                                   ,1.1 + 3 * right
                                   )#end c



                        #------------------------------------------------------------------#
                        #      95% Confidence Interval.                                    #
                        #------------------------------------------------------------------#
                        if (ed22.ci){
                           size.poly     = list()
                           size.poly$x   = c(size.poly$x,x.dbh       ,rev(x.dbh)       ,NA
                                            ,size.poly$x,x.dbh       ,rev(x.dbh)       ,NA)
                           size.poly$y   = c(size.poly$y,sta.q025[,y],rev(sta.q975[,y]),NA
                                            ,size.poly$y,ed2.q025[,y],rev(ed2.q975[,y]),NA)
                           size.poly$col = c(col.sta[2],col.ed2[2])
                        }else{
                           size.poly     = list()
                           size.poly$x   = c(size.poly$x,x.dbh       ,rev(x.dbh)       ,NA)
                           size.poly$y   = c(size.poly$y,sta.q025[,y],rev(sta.q975[,y]),NA)
                           size.poly$col = col.sta[2]
                        }#end if
                        #------------------------------------------------------------------#



                        #----- Set up the title for each plot. ----------------------------#
                        lesub = paste("Census period: ",census.desc[y],sep="")
                        #------------------------------------------------------------------#


                        #------------------------------------------------------------------#
                        #     Go on and plot stuff.                                        #
                        #------------------------------------------------------------------#
                        #----- Plotting window and grid. ----------------------------------#
                        par(mar=mar.now)
                        plot.new()
                        plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                        if (bottom) axis(side=1,at=x.dbh,labels=dbh.names)
                        if (left  ) axis(side=2)
                        box()
                        title(main=lesub)
                        if (plotgrid){
                           abline(v=x.edge,h=axTicks(2),col=grid.colour,lty="solid")
                        }#end if (plotgrid)
                        #----- Plot the taxon rate with confidence interval. --------------#
                        epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                                ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
                        lines(x=x.dbh,y=sta.expected[,y],type="o",col=col.sta[1],pch=16
                             ,lwd=2.0)
                        lines(x=x.dbh,y=ed2.expected[,y],type="o",col=col.ed2[1],pch=16
                             ,lwd=2.0)
                     }#end for
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Make the title and axis labels.                                 #
                     #---------------------------------------------------------------------#
                     letitre = paste("Size-dependent ",desc.rate,"\n",lieu,sep="")
                     ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                     lex     = desc.unit(desc="DBH class",unit=untab$cm)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the global title.                                          #
                     #---------------------------------------------------------------------#
                     gtitle( main      = letitre
                           , xlab      = lex
                           , ylab      = ley
                           , off.xlab  = 1/6
                           , line.xlab = 4.1
                           , line.ylab = 2.6
                           , cex.main  = 1.1*cex.ptsz
                           )#end gtitle
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     #---------------------------------------------------------------------#

                  }#end for
                  #------------------------------------------------------------------------#

               }#end for (i in 1:nindiv)
               #---------------------------------------------------------------------------#
            }#end for
            #==============================================================================#
            #==============================================================================#





            #==============================================================================#
            #==============================================================================#
            #  7. Plot the theme time series of the rates (by DBH class if applicable).    #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratetheme){
               this.plot  = pratetheme[[n]]
               sizetoo    = this.plot$sizetoo
               desc.rate  = this.plot$desc.rate
               unit.rate  = this.plot$unit.rate
               col.ed2    = this.plot$col.ed2
               col.sta    = this.plot$col.sta
               indiv      = this.plot$indiv
               desc.indiv = this.plot$desc.indiv
               angle      = this.plot$angle
               dens       = this.plot$density
               theme.now  = this.plot$theme
               theme.desc = this.plot$theme.desc
               plog       = this.plot$plog
               ylog       = length(grep("y",plog)) > 0

               nindiv     = length(indiv)
               nrate      = length(this.plot$ed2.rate)

               cat(" + Plotting time series of ",theme.desc,"...","\n")
               #---------------------------------------------------------------------------#
               #    Loop over all possible types of population count.                      #
               #---------------------------------------------------------------------------#
               for (i in 1:nindiv){
                  #========================================================================#
                  #========================================================================#
                  #     PLOT-LEVEL rates.                                                  #
                  #------------------------------------------------------------------------#
                  cat("  - Plot level: ",desc.indiv[i],"...","\n")
                  ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
                  sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


                  #----- Create a directory for this type of plot. ------------------------#
                  outtheme   = paste(outplot,theme.now,sep="/")
                  if (! file.exists(outtheme)) dir.create(outtheme)
                  #------------------------------------------------------------------------#


                  #----- Create path for this individual. ---------------------------------#
                  outindiv   = paste(outtheme,indiv[i],sep="/")
                  if (! file.exists(outindiv)) dir.create(outindiv)
                  #------------------------------------------------------------------------#


                  #----- Load the modelled rates. -----------------------------------------#
                  sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                  sta.expected = 100. * sta.mod[1,-1]
                  sta.q025     = 100. * sta.mod[2,-1]
                  sta.q975     = 100. * sta.mod[3,-1]
                  ed2.expected = list()
                  ed2.q025     = list()
                  ed2.q975     = list()
                  for (r in sequence(nrate)){
                     ed2.mod         = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                     ed2.expected[[r]] = 100. * ed2.mod[1,-1]
                     ed2.q025    [[r]] = 100. * ed2.mod[2,-1]
                     ed2.q975    [[r]] = 100. * ed2.mod[3,-1]
                  }#end for
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #    Find the DBH for x scale.                                           #
                  #------------------------------------------------------------------------#
                  x.years  = year4[2:n.census]
                  xlimit   = range(x.years)
                  #------------------------------------------------------------------------#


                  if (ed22.ci){
                     #---------------------------------------------------------------------#
                     #    Make the polygons.                                               #
                     #---------------------------------------------------------------------#
                     plot.poly     = list()
                     plot.poly$x   = c(x.years ,rev(x.years) )
                     plot.poly$y   = c(sta.q025,rev(sta.q975))
                     plot.poly$col = c(col.sta[2])
                     for (r in sequence(nrate)){
                        plot.poly$x   = c(plot.poly$x  ,NA
                                         ,x.years      ,rev(x.years)      )
                        plot.poly$y   = c(plot.poly$y  ,NA
                                         ,ed2.q025[[r]],rev(ed2.q975[[r]]))
                        plot.poly$col = c(plot.poly$col,col.ed2[r,2])
                     }#end for
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Find the plot limit for the y scale.                            #
                     #---------------------------------------------------------------------#
                     yuse   = c(sta.q025,sta.q975,unlist(ed2.q025),unlist(ed2.q975))
                     ylimit = pretty.xylim(u=yuse,fracexp=0.,is.log=ylog)
                     #---------------------------------------------------------------------#
                  }else{
                     #---------------------------------------------------------------------#
                     #    Make the polygons.                                               #
                     #---------------------------------------------------------------------#
                     plot.poly     = list()
                     plot.poly$x   = c(x.years ,rev(x.years) )
                     plot.poly$y   = c(sta.q025,rev(sta.q975))
                     plot.poly$col = col.sta[2]
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Find the plot limit for the y scale.                            #
                     #---------------------------------------------------------------------#
                     yuse   = c(sta.q025,sta.q975,unlist(ed2.expected))
                     ylimit = pretty.xylim(u=yuse,fracexp=0.,is.log=ylog)
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over all formats, and make the plots.                         #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Open the file or the plot window. -----------------------------#
                     fichier = file.path(outindiv
                                        ,paste(iata,"-theme-",theme.now,"-",indiv[i]
                                              ,"-",iphen.key[ph],".",outform[o],sep=""))
                     if(outform[o] == "x11"){
                        X11(width=size$width,height=size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=size$width*depth
                           ,height=size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=size$width,height=size$height
                                  ,pointsize=ptsz,paper=size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
                           ,pointsize=ptsz,paper=size$paper)
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Make the title and axis labels.                                 #
                     #---------------------------------------------------------------------#
                     letitre = paste(theme.desc," - ",lieu,sep="")
                     ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                     lex     = desc.unit(desc="Census year",unit=NULL)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Go on and plot stuff.                                           #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     layout(mat=rbind(2,1),heights=c(5,1))
                     #---------------------------------------------------------------------#


                     #----- Plot legend. --------------------------------------------------#
                     if (ed22.ci){
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0.01
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta[2],rep(0,times=nrate-1),col.ed2[,2])
                               , border  = c(col.sta[2],rep(0,times=nrate-1),col.ed2[,2])
                               , col     = c(col.sta[1],rep(0,times=nrate-1),col.ed2[,1])
                               , lwd     = c(2.0,rep(0,times=nrate-1),rep(2,0,times=nrate))
                               , pt.cex  = c(1.0,rep(0,times=nrate-1),rep(1.0,times=nrate))
                               , angle   = c(90,rep(0,times=nrate-1),angle)
                               , density = c(40,rep(0,times=nrate-1),dens)
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 1.0
                               )#end legend
                     }else{
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0.01
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta[2],rep(0,times=nrate+1))
                               , border  = c(col.sta[2],rep(0,times=nrate+1))
                               , col     = c(col.sta[1],rep(0,times=nrate-1),col.ed2[,1])
                               , lwd     = c(2.0,rep(0,times=nrate-1),rep(2,0,times=nrate))
                               , pt.cex  = c(1.0,rep(0,times=nrate-1),rep(1.0,times=nrate))
                               , angle   = c(90,rep(0,times=nrate+1))
                               , density = c(40,rep(0,times=nrate+1))
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 0.85
                               )#end legend
                     }#end if
                     #---------------------------------------------------------------------#


                     #----- Plotting window and grid. -------------------------------------#
                     par(mar=c(5,4,4,2)+0.1)
                     plot(x=x.years,y=sta.expected,xlim=xlimit,ylim=ylimit,type="n"
                         ,main=letitre,xlab=lex,ylab=ley,log=plog,cex.main=0.7)
                     if (plotgrid) grid(col=grid.colour,lty="solid")
                     #----- Plot the taxon rate with confidence interval. -----------------#
                     epolygon(x=plot.poly$x,y=plot.poly$y,col=plot.poly$col
                             ,angle=c(90,angle),density=c(40,dens),lty="solid",lwd=1.0)
                     lines(x=x.years,y=sta.expected,type="o",col=col.sta[1],pch=16,lwd=2.0)
                     for (r in 1:nrate){
                        lines(x=x.years,y=ed2.expected[[r]],type="o",col=col.ed2[r,1]
                             ,pch=16,lwd=2.0)
                     }#end for
                     #---------------------------------------------------------------------#


                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#




                  #========================================================================#
                  #========================================================================#
                  #     DBH-LEVEL rates.                                                   #
                  #------------------------------------------------------------------------#
                  if (sizetoo){
                     cat("  - DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(outsize,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Find the DBH for x scale.                                        #
                     #---------------------------------------------------------------------#
                     x.years  = year4[2:n.census]
                     xlimit   = range(x.years)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Loop over all formats, and make the plots.                      #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Open the file or the plot window. --------------------------#
                        fichier = file.path(outindiv
                                           ,paste(iata,"-theme-",theme.now,"-",indiv[i],"-"
                                                 ,iphen.key[ph],".",outform[o] ,sep=""))
                        if(outform[o] == "x11"){
                           X11(width=wide.size$width,height=wide.size$height
                              ,pointsize=ptsz)
                        }else if(outform[o] == "png"){
                           png(filename=fichier,width=wide.size$width*depth
                              ,height=wide.size$height*depth,pointsize=ptsz,res=depth)
                        }else if(outform[o] == "eps"){
                           postscript(file=fichier,width=wide.size$width
                                     ,height=wide.size$height,pointsize=ptsz
                                     ,paper=wide.size$paper)
                        }else if(outform[o] == "pdf"){
                           pdf(file=fichier,onefile=FALSE
                              ,width=wide.size$width,height=wide.size$height
                              ,pointsize=ptsz,paper=wide.size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Split the window into several smaller windows.               #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par.orig = par(no.readonly = TRUE)
                        par(oma = c(0.2,3,4,0))
                        layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                              ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#



                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        if (ed22.ci){
                           legend ( x       = "bottom"
                                  , legend  = c("Observed",rep("",times=nrate-1)
                                               ,desc.rate)
                                  , fill    = c(col.sta[2],rep(0,times=nrate-1)
                                               ,col.ed2[,2])
                                  , border  = c(col.sta[2],rep(0,times=nrate-1)
                                               ,col.ed2[,2])
                                  , col     = c(col.sta[1],rep(0,times=nrate-1)
                                               ,col.ed2[,1])
                                  , lwd     = c(2.0,rep(0,times=nrate-1)
                                                   ,rep(2,0,times=nrate))
                                  , pt.cex  = c(1.0,rep(0,times=nrate-1)
                                                   ,rep(1.0,times=nrate))
                                  , angle   = c(90,rep(0,times=nrate-1),angle)
                                  , density = c(40,rep(0,times=nrate-1),dens)
                                  , ncol    = 2
                                  , title   = "(Shaded - 95% C.I.)"
                                  , cex     = 0.85
                                  , xpd     = TRUE
                                  )#end legend
                        }else{
                           legend ( x       = "bottom"
                                  , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                                  , fill    = c(col.sta[2],rep(0,times=nrate+1))
                                  , border  = c(col.sta[2],rep(0,times=nrate+1))
                                  , col     = c(col.sta[1],rep(0,times=nrate-1)
                                               ,col.ed2[,1])
                                  , lwd     = c(2.0,rep(0,times=nrate-1)
                                                   ,rep(2,0,times=nrate))
                                  , pt.cex  = c(1.0,rep(0,times=nrate-1)
                                                   ,rep(1.0,times=nrate))
                                  , angle   = c(90,rep(0,times=nrate+1))
                                  , density = c(40,rep(0,times=nrate+1))
                                  , ncol    = 2
                                  , title   = "(Shaded - 95% C.I.)"
                                  , cex     = 0.85
                                  , xpd     = TRUE
                                  )#end legend
                        }#end if
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #    Loop over all DBH classes.                                    #
                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           left    = (d %% lo.box$ncol) == 1
                           right   = (d %% lo.box$ncol) == 0
                           top     = d <= lo.box$ncol
                           bottom  = d > (lo.box$nrow - 1) * lo.box$ncol




                           #----- Load the modelled rates. --------------------------------#
                           sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                           sta.expected = 100. * sta.mod[1,d,-1]
                           sta.q025     = 100. * sta.mod[2,d,-1]
                           sta.q975     = 100. * sta.mod[3,d,-1]
                           ed2.expected = list()
                           ed2.q025     = list()
                           ed2.q975     = list()
                           for (r in 1:nrate){
                              ed2.mod             = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                              ed2.expected[[r]]   = 100. * ed2.mod[1,d,-1]
                              ed2.q025    [[r]]   = 100. * ed2.mod[2,d,-1]
                              ed2.q975    [[r]]   = 100. * ed2.mod[3,d,-1]
                           }#end for
                           #---------------------------------------------------------------#



                           if (ed22.ci){
                              #------------------------------------------------------------#
                              #     Find the plot limit for the y scale.                   #
                              #------------------------------------------------------------#
                              yuse   = c(unlist(ed2.q025),unlist(ed2.q975)
                                        ,sta.q025,sta.q975)
                              ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #    Make the polygons.                                      #
                              #------------------------------------------------------------#
                              size.poly  = list()
                              size.poly$x   = c(x.years ,rev(x.years) )
                              size.poly$y   = c(sta.q025,rev(sta.q975))
                              size.poly$col = c(col.sta[2]  )
                              for (r in sequence(nrate)){
                                 plot.poly$x   = c(plot.poly$x  ,NA
                                                  ,x.years      ,rev(x.years)      )
                                 plot.poly$y   = c(plot.poly$y  ,NA
                                                  ,ed2.q025[[r]],rev(ed2.q975[[r]]))
                                 plot.poly$col = c(plot.poly$col,col.ed2[r,2])
                              }#end for
                              #------------------------------------------------------------#
                           }else{
                              #------------------------------------------------------------#
                              #     Find the plot limit for the y scale.                   #
                              #------------------------------------------------------------#
                              yuse   = c(ed2.expected[d,],sta.q025,sta.q975)
                              ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #    Make the polygons.                                      #
                              #------------------------------------------------------------#
                              size.poly     = list()
                              size.poly$x   = c(x.years   ,rev(x.years)     )
                              size.poly$y   = c(sta.q025  ,rev(sta.q975)    )
                              size.poly$col = c(col.sta[2])
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#


                           #----- Set up the title and axes labels. -----------------------#
                           lesub = paste("DBH class:",dbh.names[d],sep="")
                           #---------------------------------------------------------------#


                           #----- Plot the box plot. --------------------------------------#
                           par(mar=c(2,2,4,1)+0.1)
                           #----- Plotting window and grid. -------------------------------#
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           axis(side=1)
                           axis(side=2)
                           box()
                           title(main=lesub,xlab="",ylab="")
                           if (plotgrid) grid(col=grid.colour,lty="solid")
                           #----- Plot the taxon rate with confidence interval. -----------#
                           epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                                   ,angle=c(90,angle),density=c(40,dens),lty="solid"
                                   ,lwd=1.0)
                           lines(x=x.years,y=sta.expected,type="o",pch=16,lwd=2.0
                                ,col=col.sta[1])
                           for (r in sequence(nrate)){
                              lines(x=x.years,y=ed2.expected[[r]],type="o",pch=16,lwd=2.0
                                   ,col=col.ed2[r,1])
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (d in 1:n.dbh)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc,": ",lieu,sep="")
                        ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                        lex     = desc.unit(desc="Census year",unit=NULL)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Plot title.                                                  #
                        #------------------------------------------------------------------#
                        gtitle( main      = letitre
                              , xlab      = lex
                              , ylab      = ley
                              , off.xlab  = 1/6
                              , line.xlab = 4.1
                              , line.ylab = 2.6
                              , cex.main  = 1.1*cex.ptsz
                              )#end gtitle
                        #------------------------------------------------------------------#



                        #----- Close the device. ------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        #------------------------------------------------------------------#
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #==============================================================================#
            #==============================================================================#





            #==============================================================================#
            #==============================================================================#
            #  8. Plot the time series of the rates (by DBH class if applicable).          #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratets){
               this.plot  = pratets[[n]]
               sizetoo    = this.plot$sizetoo
               desc.rate  = this.plot$desc.rate
               unit.rate  = this.plot$unit.rate
               col.ed2    = this.plot$col.ed2
               col.sta    = this.plot$col.sta
               indiv      = this.plot$indiv
               desc.indiv = this.plot$desc.indiv
               colnow     = this.plot$colour
               legpos     = this.plot$legpos
               plog       = this.plot$plog
               ylog       = length(grep("y",plog)) > 0

               nindiv     = length(indiv)

               cat(" + Plotting time series of ",desc.rate,"...","\n")
               #---------------------------------------------------------------------------#
               #    Loop over all possible types of population count.                      #
               #---------------------------------------------------------------------------#
               for (i in 1:nindiv){
                  #========================================================================#
                  #========================================================================#
                  #     PLOT-LEVEL rates.                                                  #
                  #------------------------------------------------------------------------#
                  cat("  - Plot level: ",desc.indiv[i],"...","\n")
                  ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
                  sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


                  #----- Create a directory for this type of plot. ------------------------#
                  outrate   = paste(outplot,ed2.rate,sep="/")
                  if (! file.exists(outrate)) dir.create(outrate)
                  #------------------------------------------------------------------------#


                  #----- Create path for this individual. ---------------------------------#
                  outindiv   = paste(outrate,indiv[i],sep="/")
                  if (! file.exists(outindiv)) dir.create(outindiv)
                  #------------------------------------------------------------------------#


                  #----- Load the modelled rates. -----------------------------------------#
                  sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                  sta.expected = 100. * sta.mod[1,-1]
                  sta.q025     = 100. * sta.mod[2,-1]
                  sta.q975     = 100. * sta.mod[3,-1]
                  ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
                  ed2.expected = 100. * ed2.mod[1,-1]
                  ed2.q025     = 100. * ed2.mod[2,-1]
                  ed2.q975     = 100. * ed2.mod[3,-1]
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #    Find the DBH for x scale.                                           #
                  #------------------------------------------------------------------------#
                  x.years  = year4[2:n.census]
                  xlimit   = range(x.years)
                  #------------------------------------------------------------------------#


                  if (ed22.ci){
                     #---------------------------------------------------------------------#
                     #    Make the polygons.                                               #
                     #---------------------------------------------------------------------#
                     plot.poly     = list()
                     plot.poly$x   = c(x.years ,rev(x.years) ,NA,x.years ,rev(x.years) )
                     plot.poly$y   = c(sta.q025,rev(sta.q975),NA,ed2.q025,rev(ed2.q975))
                     plot.poly$col = c(col.sta[2],col.ed2[2])
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Find the plot limit for the y scale.                            #
                     #---------------------------------------------------------------------#
                     yuse   = c(sta.q025,sta.q975,ed2.q025,ed2.q975)
                     ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)
                     #---------------------------------------------------------------------#
                  }else{
                     #---------------------------------------------------------------------#
                     #    Make the polygons.                                               #
                     #---------------------------------------------------------------------#
                     plot.poly     = list()
                     plot.poly$x   = c(x.years ,rev(x.years) )
                     plot.poly$y   = c(sta.q025,rev(sta.q975))
                     plot.poly$col = col.sta[2]
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Find the plot limit for the y scale.                            #
                     #---------------------------------------------------------------------#
                     yuse   = c(sta.q025,sta.q975,ed2.expected)
                     ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over all formats, and make the plots.                         #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Open the file or the plot window. -----------------------------#
                     fichier = file.path(outindiv
                                        ,paste(iata,"-tseries-",ed2.rate,"-",indiv[i]
                                              ,"-",iphen.key[ph],".",outform[o],sep=""))
                     if(outform[o] == "x11"){
                        X11(width=size$width,height=size$height,pointsize=ptsz)
                     }else if(outform[o] == "png"){
                        png(filename=fichier,width=size$width*depth
                           ,height=size$height*depth,pointsize=ptsz,res=depth)
                     }else if(outform[o] == "eps"){
                        postscript(file=fichier,width=size$width,height=size$height
                                  ,pointsize=ptsz,paper=size$paper)
                     }else if(outform[o] == "pdf"){
                        pdf(file=fichier,onefile=FALSE
                           ,width=size$width,height=size$height,pointsize=ptsz
                           ,paper=size$paper)
                     }#end if
                     #---------------------------------------------------------------------#
                     
                     par(par.user)
                     layout(mat=rbind(2,1),heights=c(5,1))



                     #----- Plot legend. --------------------------------------------------#
                     if (ed22.ci){
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , legend  = c("Census","ED-2.2")
                               , fill    = c(col.sta[2],col.ed2[2])
                               , border  = c(col.sta[2],col.ed2[2])
                               , col     = c(col.sta[1],col.ed2[1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(-45,45)
                               , density = c( 40,40)
                               , bg      = background
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 0.85
                               , xpd     = TRUE
                               )#end legend
                     }else{
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , legend  = c("Census","ED-2.2")
                               , fill    = c(col.sta[2],         0)
                               , border  = c(col.sta[2],         0)
                               , col     = c(col.sta[1],col.ed2[1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(-45,45)
                               , density = c( 40,40)
                               , bg      = background
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 0.85
                               , xpd     = TRUE
                               )#end legend
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Make the title and axis labels.                                 #
                     #---------------------------------------------------------------------#
                     letitre = paste(desc.rate," - ",lieu,sep="")
                     ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                     lex     = desc.unit(desc="Census year",unit=NULL)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Go on and plot stuff.                                           #
                     #---------------------------------------------------------------------#
                     par(mar=c(5,4,4,2)+0.1)
                     #----- Plotting window and grid. -------------------------------------#
                     plot(x=x.years,y=sta.expected,xlim=xlimit,ylim=ylimit,type="n"
                         ,main=letitre,xlab=lex,ylab=ley,log=plog,cex.main=0.7)
                     if (plotgrid) grid(col=grid.colour,lty="solid")
                     #----- Plot the taxon rate with confidence interval. -----------------#
                     epolygon(x=plot.poly$x,y=plot.poly$y,col=plot.poly$col,angle=c(-45,45)
                             ,density=40,lty="solid",lwd=1.0)
                     lines(x=x.years,y=sta.expected,type="o",col=col.sta[1],pch=16,lwd=2.0)
                     lines(x=x.years,y=ed2.expected,type="o",col=col.ed2[1],pch=16,lwd=2.0)


                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#




                  #========================================================================#
                  #========================================================================#
                  #     DBH-LEVEL rates.                                                   #
                  #------------------------------------------------------------------------#
                  if (sizetoo){
                     cat("  - DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outrate   = paste(outsize,ed2.rate,sep="/")
                     if (! file.exists(outrate)) dir.create(outrate)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outrate,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#


                     #----- Load the modelled rates. --------------------------------------#
                     sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                     sta.expected = 100. * sta.mod[1,,-1]
                     sta.q025     = 100. * sta.mod[2,,-1]
                     sta.q975     = 100. * sta.mod[3,,-1]
                     ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
                     ed2.expected = 100. * ed2.mod[1,,-1]
                     ed2.q025     = 100. * ed2.mod[2,,-1]
                     ed2.q975     = 100. * ed2.mod[3,,-1]
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Find the DBH for x scale.                                        #
                     #---------------------------------------------------------------------#
                     x.years  = year4[2:n.census]
                     xlimit   = range(x.years)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Loop over all formats, and make the plots.                      #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Open the file or the plot window. --------------------------#
                        fichier = file.path(outindiv
                                           ,paste(iata,"-tseries-",ed2.rate,"-",indiv[i]
                                                 ,"-",iphen.key[ph],".",outform[o],sep=""))
                        if(outform[o] == "x11"){
                           X11(width=wide.size$width,height=wide.size$height
                              ,pointsize=ptsz)
                        }else if(outform[o] == "png"){
                           png(filename=fichier,width=wide.size$width*depth
                              ,height=wide.size$height*depth,pointsize=ptsz,res=depth)
                        }else if(outform[o] == "eps"){
                           postscript(file=fichier,width=wide.size$width
                                     ,height=wide.size$height,pointsize=ptsz
                                     ,paper=wide.size$paper)
                        }else if(outform[o] == "pdf"){
                           pdf(file=fichier,onefile=FALSE
                              ,width=wide.size$width,height=wide.size$height
                              ,pointsize=ptsz,paper=wide.size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Split the window into several smaller windows.               #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par.orig = par(no.readonly = TRUE)
                        par(oma = c(0.2,3,4,0))
                        layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                              ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#



                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        if (ed22.ci){
                           legend ( x       = "bottom"
                                  , inset   = 0.01
                                  , legend  = c("Census","ED-2.2")
                                  , fill    = c(col.sta[2],col.ed2[2])
                                  , border  = c(col.sta[2],col.ed2[2])
                                  , col     = c(col.sta[1],col.ed2[1])
                                  , lwd     = 2.0
                                  , pt.cex  = 1.0
                                  , angle   = c(-45,45)
                                  , density = c( 40,40)
                                  , bg      = background
                                  , ncol    = 2
                                  , title   = "(Shaded - 95% C.I.)"
                                  , cex     = 0.85
                                  )#end legend
                        }else{
                           legend ( x       = "bottom"
                                  , inset   = 0.01
                                  , legend  = c("Census","ED-2.2")
                                  , fill    = c(col.sta[2],         0)
                                  , border  = c(col.sta[2],         0)
                                  , col     = c(col.sta[1],col.ed2[1])
                                  , lwd     = 2.0
                                  , pt.cex  = 1.0
                                  , angle   = c(-45,45)
                                  , density = c( 40,40)
                                  , bg      = background
                                  , ncol    = 2
                                  , title   = "(Shaded - 95% C.I.)"
                                  , cex     = 0.85
                                  )#end legend
                        }#end if
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #    Loop over all DBH classes.                                    #
                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           left    = (d %% lo.box$ncol) == 1
                           right   = (d %% lo.box$ncol) == 0
                           top     = d <= lo.box$ncol
                           bottom  = d > (lo.box$nrow - 1) * lo.box$ncol

                           if (ed22.ci){
                              #------------------------------------------------------------#
                              #     Find the plot limit for the y scale.                   #
                              #------------------------------------------------------------#
                              yuse   = c(ed2.q025[d,],ed2.q975[d,]
                                        ,sta.q025[d,],sta.q975[d,])
                              ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #    Make the polygons.                                      #
                              #------------------------------------------------------------#
                              size.poly  = list()
                              size.poly$x   = c(x.years     ,rev(x.years)     ,NA
                                               ,x.years     ,rev(x.years)         )
                              size.poly$y   = c(sta.q025[d,],rev(sta.q975[d,]),NA
                                               ,ed2.q025[d,],rev(ed2.q975[d,])    )
                              size.poly$col = c(col.sta[2]  ,col.ed2[2]           )
                              #------------------------------------------------------------#
                           }else{
                              #------------------------------------------------------------#
                              #     Find the plot limit for the y scale.                   #
                              #------------------------------------------------------------#
                              yuse   = c(ed2.expected[d,],sta.q025[d,],sta.q975[d,])
                              ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #    Make the polygons.                                      #
                              #------------------------------------------------------------#
                              size.poly     = list()
                              size.poly$x   = c(x.years     ,rev(x.years)         )
                              size.poly$y   = c(sta.q025[d,],rev(sta.q975[d,])    )
                              size.poly$col = c(col.sta[2]  )
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#


                           #----- Set up the title and axes labels. -----------------------#
                           lesub = paste("DBH class:",dbh.names[d],sep="")
                           #---------------------------------------------------------------#


                           #----- Plot the box plot. --------------------------------------#
                           par(mar=c(2,2,4,1)+0.1)
                           #----- Plotting window and grid. -------------------------------#
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           axis(side=1)
                           axis(side=2)
                           box()
                           title(main=lesub,xlab="",ylab="")
                           if (plotgrid) grid(col=grid.colour,lty="solid")
                           #----- Plot the taxon rate with confidence interval. -----------#
                           epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                                   ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
                           lines(x=x.years,y=sta.expected[d,],type="o",pch=16,lwd=2.0
                                ,col=col.sta[1])
                           lines(x=x.years,y=ed2.expected[d,],type="o",pch=16,lwd=2.0
                                ,col=col.ed2[1])
                           #---------------------------------------------------------------#
                        }#end for (d in 1:n.dbh)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(desc.rate,": ",lieu,sep="")
                        ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                        lex     = desc.unit(desc="Census year",unit=untab$cm)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Split the plotting window.                                   #
                        #------------------------------------------------------------------#
                        gtitle( main      = letitre
                              , xlab      = lex
                              , ylab      = ley
                              , off.xlab  = 1/6
                              , line.xlab = 4.1
                              , line.ylab = 2.6
                              , cex.main  = 1.1*cex.ptsz
                              )#end gtitle
                        #------------------------------------------------------------------#



                        #----- Close the device. ------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        #------------------------------------------------------------------#
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #==============================================================================#
            #==============================================================================#
         }#end if (census.name %in% ls())
         #=================================================================================#
         #=================================================================================#










         #---------------------------------------------------------------------------------#
         #      Make the RData file name,                                                  #
         #---------------------------------------------------------------------------------#
         if (! file.exists(outdata)) dir.create(outdata)
         ed22.rdata = file.path(outdata
                               ,paste("census_",iata,"_",iphen.key[ph],"_",tfall.key[tf]
                                     ,".RData",sep=""))
         assign(x=iata,value=ed2)
         save(list=c(iata),file=ed22.rdata)
         rm(ed2)
         #---------------------------------------------------------------------------------#
      }#end for phenology
      #------------------------------------------------------------------------------------#
   }#end for tfall
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#
