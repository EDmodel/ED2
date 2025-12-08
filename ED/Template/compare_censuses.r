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


#----- Options that make path names change. -----------------------------------------------#
ibackground = 0       # Sought background colour (actual background will be transparent, 
                      #     but foreground colours will change)
                      #   0 -- white background
                      #   1 -- black background
                      #   2 -- dark grey background
oldgrowth   = FALSE   # Use old growth patches only? (FALSE uses everything)
#------------------------------------------------------------------------------------------#



#----- Main path for output. --------------------------------------------------------------#
outprefix   = ifelse(oldgrowth,"oldgrowth_comp_ibg","census_comp_ibg")
outroot     = file.path(here,paste(outprefix,sprintf("%2.2i",ibackground),sep=""))
reload      = FALSE # Should I reload data?
#------------------------------------------------------------------------------------------#



#----- Additional settings. ---------------------------------------------------------------#
yearend        = 2099
iphen.key      = c("iphen-01","iphen+02")
iphen.desc     = c("Evergreen","Drought deciduous")
tfall.key      = c("tfall111","tfall125","tfall140")[1]
tfall.desc     = c("Treefall = 1.11%/yr","Treefall = 1.25%/yr","Treefall = 1.40%/yr")[1]
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = c("pdf")               # Formats for output file.  Supported formats are:
                                        #   - "X11" - for printing on screen
                                        #   - "eps" - for postscript printing
                                        #   - "png" - for PNG printing
                                        #   - "pdf" - for PDF printing
depth          = 96                     # PNG resolution, in pixels per inch
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 18                     # Font size.
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
idbh.type      =    2         # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                              # 3 -- 0-10; 10-35; 35-55; > 55 (cm)
ed22.ci        = TRUE         # Plot confidence interval for ED?
global.ylim    = TRUE         # Global limits
n.boot         = 1000         # Number of realisations for bootstrap
#------------------------------------------------------------------------------------------#


#----- Simulation settings. ---------------------------------------------------------------#
place = list()
year.list = sprintf("%+2.2i",seq(from=-8,to=-4,by=1))
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
pop.vars      = c("n","agb","ba","acc")
pop.labels    = c("Individuals","Above Ground Biomass","Basal area","Accumulated")
pop.units     = c(untab$pcpopoyr,untab$pcagboyr,untab$pcbaoyr,untab$kgcom2oyr)
growth.vars   = c("dbh","agb","ba","acc")
growth.labels = c("DBH","Above Ground Biomass","Basal Area","Accumulated")
growth.units  = c(untab$pcdbhoyr,untab$pcagboyr,untab$pcbaoyr,untab$kgcom2oyr)
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
                   , unit.rate  = pop.units
                   , col.ed2    = c(chartreuse.fg,chartreuse.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = pop.vars
                   , desc.indiv = pop.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[2]] = list( ed2.rate   = "mort"
                   , sta.rate   = "mort"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Mortality rate"
                   , unit.rate  = pop.units
                   , col.ed2    = c(purple.fg,purple.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = pop.vars
                   , desc.indiv = pop.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[3]] = list( ed2.rate   = "ddmort"
                   , sta.rate   = "mort"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Density-dependent mort. rate"
                   , unit.rate  = pop.units
                   , col.ed2    = c(indigo.fg,indigo.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = pop.vars
                   , desc.indiv = pop.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list
pratets[[4]] = list( ed2.rate   = "dimort"
                   , sta.rate   = "mort"
                   , sizetoo    = TRUE
                   , pfttoo     = TRUE
                   , desc.rate  = "Density-independent mort. rate"
                   , unit.rate  = pop.units
                   , col.ed2    = c(blue.fg,blue.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = pop.vars
                   , desc.indiv = pop.labels
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
                     , unit.rate  = pop.units
                     , col.ed2    = c(chartreuse.fg,chartreuse.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = pop.vars
                     , desc.indiv = pop.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[2]] = list( ed2.rate   = "mort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Mortality rate"
                     , unit.rate  = pop.units
                     , col.ed2    = c(purple.fg,purple.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = pop.vars
                     , desc.indiv = pop.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[3]] = list( ed2.rate   = "ddmort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Density-dependent mort. rate"
                     , unit.rate  = pop.units
                     , col.ed2    = c(indigo.fg,indigo.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = pop.vars
                     , desc.indiv = pop.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[4]] = list( ed2.rate   = "dimort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Density-independent mort. rate"
                     , unit.rate  = pop.units
                     , col.ed2    = c(blue.fg,blue.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = pop.vars
                     , desc.indiv = pop.labels
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
pratetheme[[1]] = list( ed2.rate   = c("ddmort","dimort","mort")
                      , sta.rate   = "mort"
                      , sizetoo    = TRUE
                      , pfttoo     = TRUE
                      , desc.rate  = c("ED-2.2 Density-dependent"
                                      ,"ED-2.2 Density-independent"
                                      ,"ED-2.2 Total")
                      , unit.rate  = pop.units
                      , col.ed2    = rbind( c(sky.fg,sky.bg)
                                          , c(chartreuse.fg ,chartreuse.bg)
                                          , c(indigo.fg,indigo.bg)
                                          )#end rbind
                      , col.sta    = c(grey.fg,grey.bg)
                      , angle      = c(-30,30,60)
                      , density    = c(40,40,40)
                      , indiv      = pop.vars
                      , desc.indiv = pop.labels
                      , theme      = "mortality"
                      , theme.desc = "Mortality Rates"
                      , plog       = ""
                      )#end list
pratetheme[[2]] = list( ed2.rate   = "growth"
                      , sta.rate   = "growth"
                      , sizetoo    = TRUE
                      , pfttoo     = TRUE
                      , desc.rate  = c("ED-2.2")
                      , unit.rate  = growth.units
                      , col.ed2    = c(chartreuse.fg,chartreuse.bg)
                      , col.sta    = c(grey.mg,grey.bg)
                      , indiv      = growth.vars
                      , desc.indiv = growth.labels
                      , theme      = "productivity"
                      , theme.desc = "Growth rates"
                      , plog       = ""
                      )#end list
#---- 4. Plot expected values and confidence intervals for themes. ------------------------#
pratethbpw      = list()
pratethbpw[[1]] = list( ed2.rate   = c("mort","dimort","ddmort")
                      , sta.rate   = "mort"
                      , sizetoo    = TRUE
                      , pfttoo     = TRUE
                      , desc.rate  = c("ED-2.2 Total"
                                      ,"ED-2.2 Density-independent"
                                      ,"ED-2.2 Density-dependent")
                      , unit.rate  = pop.units
                      , col.ed2    = c(royalblue.mg,chartreuse.mg,indigo.mg)
                      , col.sta    = grey.mg
                      , indiv      = pop.vars
                      , desc.indiv = pop.labels
                      , theme      = "mortality"
                      , theme.desc = "Mortality Rates"
                      , plog       = ""
                      )#end list
pratethbpw[[2]] = list( ed2.rate   = "growth"
                      , sta.rate   = "growth"
                      , sizetoo    = TRUE
                      , pfttoo     = TRUE
                      , desc.rate  = c("ED-2.2")
                      , unit.rate  = growth.units
                      , col.ed2    = chartreuse.mg
                      , col.sta    = grey.mg
                      , indiv      = growth.vars
                      , desc.indiv = growth.labels
                      , theme      = "productivity"
                      , theme.desc = "Growth rates"
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
npratethbpw = length(pratethbpw)
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
npop.vars    = length(pop.vars  ) 
npop.vars    = length(pop.vars  ) 
ngrowth.vars = length(growth.vars) 
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size            = plotsize(proje=FALSE,paper=paper)
wide.size       = size
wide.size$width = 1.33 * size$width
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
if (! file.exists(outdata)) dir.create(outdata)
#------------------------------------------------------------------------------------------#


#----- List of places. --------------------------------------------------------------------#
nplaces = length(place)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (pl in sequence(nplaces)){
   #----- Build the list of site runs. ----------------------------------------------------#
   iata      = place[[pl]]$iata
   idx       = match(iata,poilist$iata)
   longname  = poilist$longname[idx]
   config    = expand.grid(place[[pl]]$config,stringsAsFactors=FALSE)
   simul.all = paste("t",iata,"_",apply(X=config,MARGIN=1,FUN=paste,collapse="_")
                    ,sep="")
   #---------------------------------------------------------------------------------------#




   #----- Find the census observations for this particular site. --------------------------#
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
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     We only run this part of the code if there are observations to compare with       #
   # the model.                                                                            #
   #---------------------------------------------------------------------------------------#
   if (census.name %in% ls()){

      #------------------------------------------------------------------------------------#
      #     Obs.rdata is similar to the input census data, but with extra variables to     #
      # control plots.                                                                     #
      #------------------------------------------------------------------------------------#
      cat (" + Site: ",longname,"\n")
      if (oldgrowth){
         ylim.rdata = file.path(outdata
                               ,paste("oldgrowth_",iata,"_","limits",".RData",sep=""))
      }else{
         ylim.rdata = file.path(outdata
                               ,paste("census_",iata,"_","limits",".RData",sep=""))
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Load the census data, from the monthly means.                                  #
      #------------------------------------------------------------------------------------#
      sta         = get(census.name)
      n.census    = length(sta$when)
      n.dbh       = length(sta$dbh.breaks)-1
      
      x.edge      = c(10,      sta$dbh.breaks[-c(1,n.dbh+1)]
                        , 2. * sta$dbh.breaks[n.dbh] - sta$dbh.breaks[n.dbh-1] )
      x.dbh       = 0.5 * ( x.edge[-1] + x.edge[-(n.dbh+1)] )
      xlimit      = pretty.xylim(u=x.edge,fracexp=0.0,is.log=FALSE)
      dbh.names   = dimnames(sta$mort.size$n$expected)[[2]]
      year4       = numyears(sta$when)
      when4       = sta$when
      when.label  = paste(month.abb[nummonths(when4)],year4,sep="-")
      when.limit  = c(chron(paste(1,1,min(year4)  ,sep="/"))
                     ,chron(paste(1,1,max(year4)+1,sep="/")))
      whenmid4    = chron(0.5 * ( as.numeric( sta$when[       -1] )
                                + as.numeric( sta$when[-n.census] ) ) )

      delta       = diff(as.numeric(sta$when))
      mid.label   = c(chron(sta$when[1]-delta[1]),sta$when)
      month.label = month.abb[nummonths(mid.label)]
      year.label  = numyears(mid.label)
      
      summ.label  = paste(month.label,year.label,sep="/")
      n.label     = length(summ.label)

      census.desc = paste(summ.label[-n.label],summ.label[-1],sep=" - ")
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #      Loop over all months to grab all the census data.                             #
      #------------------------------------------------------------------------------------#
      census.idx   = NULL
      for (y in 2:n.census){
         #----- Find the first and last time to be averaged for this census. --------------#
         ts.montha  = ( nummonths(sta$when[y-1]) %% 12 )
         ts.yeara   = numyears (sta$when[y-1])
         ts.monthz  = ( ( (nummonths(sta$when[y]) - 1) %% 12 )
                      + 12 * as.integer(nummonths(sta$when[y]) == 1) )
         ts.yearz   = numyears (sta$when[y]) - as.integer(ts.monthz == 12)
         n.inter    = (ts.yearz-ts.yeara-1)*12 + ts.monthz + (12 - ts.montha + 1)
         #---------------------------------------------------------------------------------#
         census.idx = c(census.idx,rep(y,times=n.inter))
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the average rates for the census period using log-normal.                 #
      #------------------------------------------------------------------------------------#
      mypfts = sort(match(unique(names(sta$classes)),pft$name))
      npfts  = length(mypfts)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Loop over all the possible rates and initialise the limits.                   #
      #------------------------------------------------------------------------------------#
      ylt = list()
      for (r in sequence(npratets)){
         #---------------------------------------------------------------------------------#
         #     Load the rate information.                                                  #
         #---------------------------------------------------------------------------------#
         this.rate = pratets[[ r]]
         ed2.rate  = this.rate$ed2.rate
         sta.rate  = this.rate$sta.rate
         sizetoo   = this.rate$sizetoo
         desc.rate = this.rate$desc.rate
         indiv     = this.rate$indiv
         cat(" - Initialising limits for ",desc.rate," tables...","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Load plot-level and size-level data.                                       #
         #---------------------------------------------------------------------------------#
         sta.plot   = paste(sta.rate,"plot",sep=".")
         sta.size   = paste(sta.rate,"size",sep=".")
         ed2.plot   = paste(ed2.rate,"plot",sep=".")
         ed2.size   = paste(ed2.rate,"size",sep=".")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Find how many individuals to retrieve.                                    #
         #---------------------------------------------------------------------------------#
         nindiv = length(indiv)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop over the different types of individuals.                               #
         #---------------------------------------------------------------------------------#
         ylt[[ed2.plot]] = list()
         ylt[[ed2.size]] = list()
         for (v in 1:nindiv){
            #----- Set up the individuals. ------------------------------------------------#
            vn   = indiv[v]
            yy   = 2:n.census
            #------------------------------------------------------------------------------#


            #----- Grab minima and maxima from quantiles. ---------------------------------#
            sta.now         = sta[[sta.plot]][[vn]]
            ylim.min.global = min(sta.now$global[2,],na.rm=TRUE)
            ylim.max.global = max(sta.now$global[3,],na.rm=TRUE)
            ylim.min.taxon  = apply(sta.now$q025,MARGIN=1,FUN=min,na.rm=TRUE)
            ylim.max.taxon  = apply(sta.now$q975,MARGIN=1,FUN=max,na.rm=TRUE)
            #------------------------------------------------------------------------------#



            #----- Initialise the plot-level structure. -----------------------------------#
            ylt[[ed2.plot]][[vn]] = list( global = list ( min = ylim.min.global
                                                        , max = ylim.max.global
                                                        )#end list
                                        , taxon  = list ( min = ylim.min.taxon
                                                        , max = ylim.max.taxon
                                                        )#end list
                                        )#end list
            rm(sta.now)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #    Now the size-dependent variables, if this is a size-dependent rate.       #
            #------------------------------------------------------------------------------#
            if (sizetoo){


               #----- Grab minima and maxima from quantiles. ------------------------------#
               sta.now         = sta[[sta.size]][[vn]]
               ylim.min.global = apply(sta.now$global[2,,],1,min,na.rm=TRUE)
               ylim.max.global = apply(sta.now$global[3,,],1,max,na.rm=TRUE)
               ylim.min.taxon  = apply(sta.now$q025,c(1,2),FUN=min,na.rm=TRUE)
               ylim.max.taxon  = apply(sta.now$q975,c(1,2),FUN=max,na.rm=TRUE)
               #---------------------------------------------------------------------------#


               #----- Initialise the size-level structure. --------------------------------#
               ylt[[ed2.size]][[vn]] = list( global = list ( min = ylim.min.global
                                                           , max = ylim.max.global
                                                           )#end list
                                           , taxon  = list ( min   = ylim.min.taxon
                                                           , max   = ylim.max.taxon
                                                           )#end list
                                           )#end list
               rm(sta.now)
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Treefall loop.                                                                 #
      #------------------------------------------------------------------------------------#
      for (tf in sequence(n.tfall)){
         #----- Set up the path for this treefall. ----------------------------------------#
         outfall = file.path(outroot,tfall.key[tf])
         if (! file.exists(outfall)) dir.create(outfall)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Phenology loop.                                                             #
         #---------------------------------------------------------------------------------#
         for (ph in sequence(n.iphen)){

            cat ("  - Group: ",tfall.desc[tf]," - ",iphen.desc[ph],"\n")



            #----- Use only this phenology and treefall combination. ----------------------#
            aux.iphen.key=sub(pattern="\\+",replacement=".",x=iphen.key)
            aux.tfall.key=sub(pattern="\\+",replacement=".",x=tfall.key)
            aux.simul.all=sub(pattern="\\+",replacement=".",x=simul.all)
            keep    = ( regexpr(pattern=aux.iphen.key[ph],text=aux.simul.all) > 0
                      & regexpr(pattern=aux.tfall.key[tf],text=aux.simul.all) > 0 )
            simul   = simul.all[keep]
            n.simul = length(simul)
            #------------------------------------------------------------------------------#




            #------ Find the output file. -------------------------------------------------#
            if (! file.exists(outdata)) dir.create(outdata)
            if (oldgrowth){
               ed22.rdata = file.path(outdata,paste("oldgrowth_",iata,"_",iphen.key[ph],"_"
                                                   ,tfall.key[tf],".RData",sep=""))
            }else{
               ed22.rdata = file.path(outdata,paste("census_",iata,"_",iphen.key[ph],"_"
                                                   ,tfall.key[tf],".RData",sep=""))
            }#end if
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
               if (oldgrowth){
                  this.rdata = file.path(outmain,"rdata_census"
                                        ,paste("oldgrowth_",simul[s],".RData",sep=""))
               }else{
                  this.rdata = file.path(outmain,"rdata_census"
                                        ,paste("census_",simul[s],".RData",sep=""))
               }#end if
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
                  for (v in 1:npop.vars){
                     ts.recr.size  [[pop.vars[v]]] = array( NA, dim = dim.size)
                  }#end for
                  #----- Mortality. -------------------------------------------------------#
                  ts.mort.size     = list()
                  ts.ddmort.size   = list()
                  ts.dimort.size   = list()
                  for (v in 1:npop.vars){
                     ts.mort.size  [[pop.vars[v]]] = array( NA, dim = dim.size)
                     ts.ddmort.size[[pop.vars[v]]] = array( NA, dim = dim.size)
                     ts.dimort.size[[pop.vars[v]]] = array( NA, dim = dim.size)
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
               for (v in 1:npop.vars){
                  var.now = pop.vars[v]
                  ts.recr.size  [[var.now]][,,,,s] = ed2$tseries$recr  [[var.now]]
               }#end for
               for (v in 1:npop.vars){
                  var.now = pop.vars[v]
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

            #------------------------------------------------------------------------------#
            #     Loop over all rates.                                                     #
            #------------------------------------------------------------------------------#
            recr.rates = c("recr")
            n.rates    = length(recr.rates)
            for (r in 1:n.rates){
               ts.this.size = get(paste("ts",recr.rates[r],"size",sep="."))
               ms.this.plot = list()
               ms.this.size = list()

               #---------------------------------------------------------------------------#
               #     Loop over all rates.                                                  #
               #---------------------------------------------------------------------------#
               for (v in 1:npop.vars){
                  v.now  = pop.vars[v]
                  ms.now = list()


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
                        ts.plot.now       = c(ts.this.size[[v.now]][p,n.dbh+1,i.sel,,])
                        ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                        if (any(is.finite(ts.plot.now))){
                           boot.now = shhh(fun=boot,data=ts.plot.now,statistic=mean.fun
                                          ,R=n.boot)
                           ci.now   = try(shhh(fun=boot.ci,boot.out=boot.now
                                              ,conf=0.95,type="perc")
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
                              boot.now = shhh(fun=boot,data=ts.size.now,statistic=mean.fun
                                             ,R=n.boot)
                              ci.now   = try(shhh(fun=boot.ci,boot.out=boot.now
                                                 ,conf=0.95,type="perc")
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
               }#end for
               #---------------------------------------------------------------------------#
               dummy = assign(paste("ms",recr.rates[r],"plot",sep="."), ms.this.plot)
               dummy = assign(paste("ms",recr.rates[r],"size",sep="."), ms.this.size)
            }#end for
            #------------------------------------------------------------------------------#


            #----- Mortality rates. -------------------------------------------------------#
            cat("     * Mortality...","\n")

            mort.rates = c("mort","ddmort","dimort")
            n.rates    = length(mort.rates)

            #---------------------------------------------------------------------------#
            #     Loop over all rates.                                                  #
            #---------------------------------------------------------------------------#
            for (r in 1:n.rates){
               ts.this.size = get(paste("ts",mort.rates[r],"size",sep="."))
               ms.this.plot = list()
               ms.this.size = list()
               for (v in 1:npop.vars){
                  v.now  = pop.vars[v]
                  ms.now = list()


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
                        ts.plot.now       = c(ts.this.size[[v.now]][p,n.dbh+1,i.sel,,])
                        ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                        if (any(is.finite(ts.plot.now))){
                           boot.now = shhh(fun=boot,data=ts.plot.now,statistic=mean.fun
                                          ,R=n.boot)
                           ci.now   = try(shhh(fun=boot.ci,boot.out=boot.now
                                              ,conf=0.95,type="perc")
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
                              boot.now = shhh(fun=boot,data=ts.size.now,statistic=mean.fun
                                             ,R=n.boot)
                              ci.now   = try(shhh(fun=boot.ci,boot.out=boot.now
                                                 ,conf=0.95,type="perc")
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
               }#end for
               dummy = assign(paste("ms",mort.rates[r],"plot",sep="."), ms.this.plot)
               dummy = assign(paste("ms",mort.rates[r],"size",sep="."), ms.this.size)
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#

            
            #----- Growth rates. ----------------------------------------------------------#
            cat("     * Growth...","\n")

            #---------------------------------------------------------------------------#
            #     Loop over all rates.                                                  #
            #---------------------------------------------------------------------------#
            growth.rates = c("growth")
            n.rates    = length(growth.rates)
            for (r in 1:n.rates){
               ts.this.size = get(paste("ts",growth.rates[r],"size",sep="."))
               ms.this.plot = list()
               ms.this.size = list()

               for (v in 1:ngrowth.vars){
                  v.now  = growth.vars[v]
                  ms.now = list()


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
                        ts.plot.now       = c(ts.this.size[[v.now]][p,n.dbh+1,i.sel,,])
                        ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                        if (any(is.finite(ts.plot.now))){
                           boot.now = shhh(fun=boot,data=ts.plot.now,statistic=mean.fun
                                          ,R=n.boot)
                           ci.now   = try(shhh(fun=boot.ci,boot.out=boot.now
                                              ,conf=0.95,type="perc")
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
                              boot.now = shhh(fun=boot,data=ts.size.now,statistic=mean.fun
                                             ,R=n.boot)
                              ci.now   = try(shhh(fun=boot.ci,boot.out=boot.now
                                                 ,conf=0.95,type="perc")
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
               }#end for
               dummy = assign(paste("ms",growth.rates[r],"plot",sep="."), ms.this.plot)
               dummy = assign(paste("ms",growth.rates[r],"size",sep="."), ms.this.size)
               #---------------------------------------------------------------------------#
            }#end for
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
            #      Loop over all the possible rates.                                       #
            #------------------------------------------------------------------------------#
            for (r in sequence(npratets)){
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
                  #     Update limits.                                                     #
                  #------------------------------------------------------------------------#
                  ed2.now = ed2[[ed2.plot]][[vn]]
                  ylt.now = ylt[[ed2.plot]][[vn]]
                  ylim.min.global = min  (x=ed2.now$global[2,],na.rm=TRUE)
                  ylim.max.global = max  (x=ed2.now$global[3,],na.rm=TRUE)
                  ylim.min.taxon  = apply(X=ed2.now$q025,MARGIN=1,FUN=min,na.rm=TRUE)
                  ylim.max.taxon  = apply(X=ed2.now$q975,MARGIN=1,FUN=max,na.rm=TRUE)

                  ylt[[ed2.plot]][[vn]]$global$min = min ( ylt.now$global$min
                                                         , ylim.min.global
                                                         , na.rm=TRUE
                                                         )#end min
                  ylt[[ed2.plot]][[vn]]$global$max = max ( ylt.now$global$max
                                                         , ylim.max.global
                                                         , na.rm=TRUE
                                                         )#end max
                  ylt[[ed2.plot]][[vn]]$taxon$min  = pmin( ylt.now$taxon$min
                                                         , ylim.min.taxon
                                                         , na.rm = TRUE
                                                         )#end pmin
                  ylt[[ed2.plot]][[vn]]$taxon$max  = pmax( ylt.now$taxon$max
                                                         , ylim.max.taxon
                                                         , na.rm = TRUE
                                                         )#end pmax
                  rm(ed2.now,ylt.now)
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





                     #---------------------------------------------------------------------#
                     #    Update limits.                                                   #
                     #---------------------------------------------------------------------#
                     ed2.now = ed2[[ed2.size]][[vn]]
                     ylt.now = ylt[[ed2.size]][[vn]]

                     ylim.min.global = apply(X=ed2.now$global[2,,],1     ,min,na.rm=TRUE)
                     ylim.max.global = apply(X=ed2.now$global[3,,],1     ,max,na.rm=TRUE)
                     ylim.min.taxon  = apply(X=ed2.now$q025       ,c(1,2),min,na.rm=TRUE)
                     ylim.max.taxon  = apply(X=ed2.now$q975       ,c(1,2),max,na.rm=TRUE)

                     ylt[[ed2.size]][[vn]]$global$min = pmin( ylt.now$global$min
                                                            , ylim.min.global
                                                            , na.rm=TRUE
                                                            )#end min
                     ylt[[ed2.size]][[vn]]$global$max = pmax( ylt.now$global$max
                                                            , ylim.max.global
                                                            , na.rm=TRUE
                                                            )#end max
                     ylt[[ed2.size]][[vn]]$taxon$min  = pmin( ylt.now$taxon$min
                                                            , ylim.min.taxon
                                                            , na.rm = TRUE
                                                            )#end pmin
                     ylt[[ed2.size]][[vn]]$taxon$max  = pmax( ylt.now$taxon$max
                                                            , ylim.max.taxon
                                                            , na.rm = TRUE
                                                            )#end pmax
                     rm(ed2.now,ylt.now)
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Make the RData file name,                                               #
            #------------------------------------------------------------------------------#
            cat (" + Saving simulation to ",basename(ed22.rdata),"...","\n")
            assign(x=iata,value=ed2)
            save(list=c(iata),file=ed22.rdata)
            rm(ed2)
            #------------------------------------------------------------------------------#
         }#end for phenology
         #---------------------------------------------------------------------------------#
      }#end for tfall
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Make the RData file name,                                                     #
      #------------------------------------------------------------------------------------#
      cat (" + Saving limits to ",basename(ylim.rdata),"...","\n")
      assign(x=iata,value=ylt)
      save(list=c(iata),file=ylim.rdata)
      rm(ylt)
      #------------------------------------------------------------------------------------#
   }#end if (census.name %in% ls())
   #=======================================================================================#
   #=======================================================================================#
}#end for places
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#     The other big place loop starts here...                                              #
#------------------------------------------------------------------------------------------#
for (pl in sequence(nplaces)){
   #----- Build the list of site runs. ----------------------------------------------------#
   iata      = place[[pl]]$iata
   idx       = match(iata,poilist$iata)
   longname  = poilist$longname[idx]
   config    = expand.grid(place[[pl]]$config,stringsAsFactors=FALSE)
   simul.all = paste("t",iata,"_",apply(X=config,MARGIN=1,FUN=paste,collapse="_")
                    ,sep="")
   #---------------------------------------------------------------------------------------#




   #----- Find the census observations for this particular site. --------------------------#
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
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     We only run this part of the code if there are observations to compare with the   #
   # model.                                                                                #
   #---------------------------------------------------------------------------------------#
   if (census.name %in% ls()){
      #------------------------------------------------------------------------------------#
      #     Load the census data, from the monthly means.                                  #
      #------------------------------------------------------------------------------------#
      cat (" + Site: ",longname,"\n")
      sta         = get(census.name)
      n.census    = length(sta$when)
      n.dbh       = length(sta$dbh.breaks)-1
      
      x.edge      = c(10,      sta$dbh.breaks[-c(1,n.dbh+1)]
                        , 2. * sta$dbh.breaks[n.dbh] - sta$dbh.breaks[n.dbh-1] )
      x.dbh       = 0.5 * ( x.edge[-1] + x.edge[-(n.dbh+1)] )
      xlimit      = pretty.xylim(u=x.edge,fracexp=0.0,is.log=FALSE)
      dbh.names   = dimnames(sta$mort.size$n$expected)[[2]]
      year4       = numyears(sta$when)
      when4       = sta$when
      when.label  = paste(month.abb[nummonths(when4)],year4,sep="-")
      when.limit  = c(chron(paste(1,1,min(year4)  ,sep="/"))
                     ,chron(paste(1,1,max(year4)+1,sep="/")))
      whenmid4    = chron(0.5 * ( as.numeric( sta$when[       -1] )
                                + as.numeric( sta$when[-n.census] ) ) )

      delta       = diff(as.numeric(sta$when))
      mid.label   = c(chron(sta$when[1]-delta[1]),sta$when)
      month.label = month.abb[nummonths(mid.label)]
      year.label  = numyears(mid.label)
      
      summ.label  = paste(month.label,year.label,sep="/")
      n.label     = length(summ.label)

      census.desc = paste(summ.label[-n.label],summ.label[-1],sep=" - ")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Load limits for this place.                                                   #
      #------------------------------------------------------------------------------------#
      ylim.rdata = file.path(outdata,paste("census_",iata,"_","limits",".RData",sep=""))
      load(ylim.rdata)
      ylt = get(iata)
      rm(iata)
      iata      = place[[pl]]$iata
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Treefall loop.                                                                 #
      #------------------------------------------------------------------------------------#
      for (tf in sequence(n.tfall)){
         #----- Set up the path for this treefall. ----------------------------------------#
         outfall = file.path(outroot,tfall.key[tf])
         if (! file.exists(outfall)) dir.create(outfall)
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #     Make the directories.                                                       #
         #---------------------------------------------------------------------------------#
         outplot   = file.path(outfall,"census_plot")
         outsize   = file.path(outfall,"census_size")
         pftplot   = file.path(outfall,"pft_plot"   )
         pftsize   = file.path(outfall,"pft_size"   )
         if (! file.exists(outplot)) dir.create(outplot)
         if (! file.exists(outsize)) dir.create(outsize)
         if (! file.exists(pftplot)) dir.create(pftplot)
         if (! file.exists(pftsize)) dir.create(pftsize)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Phenology loop.                                                             #
         #---------------------------------------------------------------------------------#
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



            #------------------------------------------------------------------------------#
            #      Load current observation.                                               #
            #------------------------------------------------------------------------------#
            ed22.rdata = file.path(outdata
                                  ,paste("census_",iata,"_",iphen.key[ph],"_",tfall.key[tf]
                                        ,".RData",sep=""))
            load(ed22.rdata)
            ed2 = get(iata)
            rm(iata)
            iata      = place[[pl]]$iata
            #------------------------------------------------------------------------------#






            #==============================================================================#
            #==============================================================================#
            #  6.  Plot rates as function of size.                                         #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratesize){
               this.plot  = pratesize[[n]]
               pfttoo     = this.plot$pfttoo
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


               cat(" + Plotting size-dependent ",desc.rate,"...","\n")

               #---------------------------------------------------------------------------#
               #    Loop over all possible types of population count.                      #
               #---------------------------------------------------------------------------#
               for (i in 1:nindiv){
                  cat("  - Size level: ",desc.indiv[i],"...","\n")


                  #----- Build the rate name. ---------------------------------------------#
                  ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                  sta.rate = paste(this.plot$sta.rate,"size",sep=".")
                  #------------------------------------------------------------------------#



                  #----- Create a directory for this type of plot. ------------------------#
                  outrate   = paste(outsize,ed2.rate,sep="/")
                  if (! file.exists(outrate)) dir.create(outrate)
                  if (pfttoo){
                     pftrate   = paste(pftsize,ed2.rate,sep="/")
                     if (! file.exists(pftrate)) dir.create(pftrate)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Create path for this individual. ---------------------------------#
                  outindiv   = paste(outrate,indiv[i],sep="/")
                  if (! file.exists(outindiv)) dir.create(outindiv)
                  #------------------------------------------------------------------------#


                  #----- Load the modelled rates. -----------------------------------------#
                  mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                  sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                  sta.expected = mult * sta.mod[1,,]
                  sta.q025     = mult * sta.mod[2,,]
                  sta.q975     = mult * sta.mod[3,,]
                  ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
                  ed2.expected = mult * ed2.mod[1,,]
                  ed2.q025     = mult * ed2.mod[2,,]
                  ed2.q975     = mult * ed2.mod[3,,]
                  ylt.mod      = ylt[[ed2.rate]][[indiv[i]]]$global
                  ylt.min      = mult * ylt.mod$min
                  ylt.max      = mult * ylt.mod$max
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
                        yuse   = c(sta.q025,sta.q975,sta.expected
                                  ,ed2.q025,ed2.q975,ed2.expected)
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
                        yuse   = c(sta.q025,sta.q975,sta.expected,ed2.expected)
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
                     #     Overwrite limits in case the user wants global limits.          #
                     #---------------------------------------------------------------------#
                     if (global.ylim){
                        ylimit = pretty.xylim(u=c(ylt.min,ylt.max),fracexp=0.0,is.log=ylog)
                     }#end if
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Loop over all years.                                           #
                     #---------------------------------------------------------------------#
                     for (y in 2:n.census){
                        k       = y - 1
                        left    = lo.box$left  [k ]
                        bottom  = lo.box$bottom[k ]
                        mar.now = lo.box$mar   [k,]



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
                        if (left  ) axis(side=2,las=1)
                        box()
                        title(main=lesub)
                        if (plotgrid){
                           abline(v=x.edge,h=axTicks(2),col=grid.colour
                                 ,lwd=0.75,lty="dotted")
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
                     letitre = paste("Size-dependent ",desc.rate," - ",longname
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
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
                           , cex.main  = 0.9*cex.ptsz
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






                  #------------------------------------------------------------------------#
                  #     Plot the size structure by PFT.                                    #
                  #------------------------------------------------------------------------#
                  if (pfttoo){
                     cat("   - Size level by PFT: ",desc.indiv[i]," ...","\n")

                     #----- Build the rate name. ------------------------------------------#
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")
                     #---------------------------------------------------------------------#



                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(pftrate,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#


                     for (p in sequence(npfts)){
                        pft.label = pft$name[mypfts[p]]
                        pft.key   = paste("pft",sprintf("%2.2i",mypfts[p]),sep="")


                        #----- Load the modelled rates. -----------------------------------#
                        mult         = 100 - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                        sta.mod      = sta[[sta.rate]][[indiv[i]]]
                        sta.expected = mult * sta.mod$expected [p,,]
                        sta.q025     = mult * sta.mod$q025     [p,,]
                        sta.q975     = mult * sta.mod$q975     [p,,]
                        ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]
                        ed2.expected = mult * ed2.mod$expected [p,,]
                        ed2.q025     = mult * ed2.mod$q025     [p,,]
                        ed2.q975     = mult * ed2.mod$q975     [p,,]
                        ylt.mod      = ylt[[ed2.rate]][[indiv[i]]]
                        ylt.min      = mult * ylt.mod$taxon$min[p, ]
                        ylt.max      = mult * ylt.mod$taxon$min[p, ]
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Define a nice configuration for the multiple panels.        #
                        #------------------------------------------------------------------#
                        lo.box = pretty.box(n=n.census-1,horizontal=TRUE)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Loop over all formats.                                      #
                        #------------------------------------------------------------------#
                        for (o in 1:nout){
                           #----- Open the file or the plot window. -----------------------#
                           fichier = file.path(outindiv
                                              ,paste(iata,"-yrsize-",ed2.rate,"-",indiv[i]
                                                    ,"-",pft.key,"-",iphen.key[ph]
                                                    ,".",outform[o],sep=""))
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
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Split the window into several smaller windows.            #
                           #---------------------------------------------------------------#
                           par(par.user)
                           par.orig = par(no.readonly = TRUE)
                           par(oma = c(0.2,3,4,0))
                           layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                                 ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                                 )#end layout
                           #---------------------------------------------------------------#


                           #---------------------------------------------------------------#
                           #     Find the plot limit for the y scale.                      #
                           #---------------------------------------------------------------#
                           if (ed22.ci){
                              yuse   = c(sta.q025,sta.q975,sta.expected
                                        ,ed2.q025,ed2.q975,ed2.expected)
                              ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)



                              #----- Plot legend. -----------------------------------------#
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
                              #------------------------------------------------------------#
                           }else{
                              yuse   = c(sta.q025,sta.q975,sta.expected,ed2.expected)
                              ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)



                              #----- Plot legend. -----------------------------------------#
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
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #      Overwrite if global.ylim is TRUE.                        #
                           #---------------------------------------------------------------#
                           if (global.ylim){
                              ylimit = pretty.xylim(u=c(ylt.min,ylt.max),fracexp=0.0
                                                   ,is.log=ylog)
                           }#end if
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #      Loop over all years.                                     #
                           #---------------------------------------------------------------#
                           for (y in 2:n.census){
                              k       = y - 1
                              left    = lo.box$left  [k ]
                              bottom  = lo.box$bottom[k ]
                              mar.now = lo.box$mar   [k,]



                              #------------------------------------------------------------#
                              #      95% Confidence Interval.                              #
                              #------------------------------------------------------------#
                              if (ed22.ci){
                                 size.poly     = list()
                                 size.poly$x   = c(size.poly$x,x.dbh       
                                                  ,rev(x.dbh)       ,NA
                                                  ,size.poly$x,x.dbh       
                                                  ,rev(x.dbh)       ,NA)
                                 size.poly$y   = c(size.poly$y,sta.q025[,y]
                                                  ,rev(sta.q975[,y]),NA
                                                  ,size.poly$y,ed2.q025[,y]
                                                  ,rev(ed2.q975[,y]),NA)
                                 size.poly$col = c(col.sta[2],col.ed2[2])
                              }else{
                                 size.poly     = list()
                                 size.poly$x   = c(size.poly$x,x.dbh       
                                                  ,rev(x.dbh)       ,NA)
                                 size.poly$y   = c(size.poly$y,sta.q025[,y]
                                                  ,rev(sta.q975[,y]),NA)
                                 size.poly$col = col.sta[2]
                              }#end if
                              #------------------------------------------------------------#



                              #----- Set up the title for each plot. ----------------------#
                              lesub = paste("Census period: ",census.desc[y],sep="")
                              #------------------------------------------------------------#


                              #------------------------------------------------------------#
                              #     Go on and plot stuff.                                  #
                              #------------------------------------------------------------#
                              #----- Plotting window and grid. ----------------------------#
                              par(mar=mar.now)
                              plot.new()
                              plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                              if (bottom) axis(side=1,at=x.dbh,labels=dbh.names)
                              if (left  ) axis(side=2,las=1)
                              box()
                              title(main=lesub)
                              if (plotgrid){
                                 abline(v=x.edge,h=axTicks(2),col=grid.colour
                                       ,lwd=0.75,lty="dotted")
                              }#end if (plotgrid)
                              #----- Plot the taxon rate with confidence interval. --------#
                              epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                                      ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
                              lines(x=x.dbh,y=sta.expected[,y],type="o",col=col.sta[1]
                                   ,pch=16,lwd=2.0)
                              lines(x=x.dbh,y=ed2.expected[,y],type="o",col=col.ed2[1]
                                   ,pch=16,lwd=2.0)
                           }#end for
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Make the title and axis labels.                           #
                           #---------------------------------------------------------------#
                           letitre = paste("Size-dependent ",desc.rate," - ",longname
                                          ," - ",pft.label
                                          ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                           ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                           lex     = desc.unit(desc="DBH class",unit=untab$cm)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Plot the global title.                                    #
                           #---------------------------------------------------------------#
                           gtitle( main      = letitre
                                 , xlab      = lex
                                 , ylab      = ley
                                 , off.xlab  = 1/6
                                 , line.xlab = 4.1
                                 , line.ylab = 2.6
                                 , cex.main  = 0.9*cex.ptsz
                                 )#end gtitle
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           #---------------------------------------------------------------#

                        }#end for (o in 1:nout)
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts))
                  }#end if (pfttoo)
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
               pfttoo     = this.plot$pfttoo
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
               col.ed2    = matrix(col.ed2,nrow=nrate)

               cat(" + Plotting theme time series of ",theme.desc,"...","\n")
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
                  mult   = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                  mdot   = list()
                  mci    = list()
                  xlimit = when.limit
                  ylimit = NULL
                  ny     = n.census
                  for (r in sequence(nrate+1)){
                     if (r == nrate+1){
                        mod.now      = sta[[sta.rate       ]][[indiv[i]]]$global
                        cimult       = mult
                        col.exp      = col.sta[1]
                        col.ci       = col.sta[2]
                        angle.ci     = 90
                        dens.ci      = 40
                     }else{
                        mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                        cimult       = ifelse(ed22.ci,mult,NA)
                        col.exp      = col.ed2[r,1]
                        col.ci       = col.ed2[r,2]
                        angle.ci     = angle[r]
                        dens.ci      = dens [r]
                     }#end if

                     #----- Store polygons. -----------------------------------------------#
                     mdot[[r]] = list( x      = whenmid4
                                     , y      = mult*mod.now[1,-1]
                                     , col    = col.exp
                                     , pch    = 16
                                     , type   = "o"
                                     )#end list
                     mci [[r]] = list( x = c(whenmid4,rev(whenmid4))
                                     , y = cimult * c(mod.now[2,-1],rev(mod.now[3,-1]))
                                     , col   = col.ci
                                     , angle = angle.ci
                                     , dens  = dens.ci
                                     , lty   = "solid"
                                     , lwd   = 1.0
                                     )#end list
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Update y range.                                                 #
                     #---------------------------------------------------------------------#
                     if (global.ylim && r != (nrate + 1)){
                        ylt.min = mult * ylt[[ed2.rate[r]]][[indiv[i]]]$global$min
                        ylt.max = mult * ylt[[ed2.rate[r]]][[indiv[i]]]$global$max
                        ylimit  = range(c(ylimit,ylt.min,ylt.max))
                     }else{
                        ylimit  = range(c(ylimit,mdot[[r]]$y,mci[[r]]$y),na.rm=TRUE)
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for (r in sequence(nrate+1))
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
                     letitre = paste(theme.desc," - ",longname
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                     ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                     lex     = "" # desc.unit(desc="Census year",unit=NULL)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Go on and plot stuff.                                           #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     layout(mat=rbind(2,1),heights=c(5,1))
                     #---------------------------------------------------------------------#


                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.1,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                            , fill    = c(col.sta[2],rep("transparent",times=nrate-1)
                                         ,col.ed2[,2])
                            , border  = c(col.sta[2],rep("transparent",times=nrate-1)
                                         ,col.ed2[,2])
                            , col     = c(col.sta[1],rep("transparent",times=nrate-1)
                                         ,col.ed2[,1])
                            , lwd     = 2.0
                            , pt.cex  = 1.0
                            , angle   = c(90,rep(0,times=nrate-1),angle)
                            , density = c(40,rep(0,times=nrate-1),dens)
                            , ncol    = 2
                            , title   = "(Shaded - 95% C.I.)"
                            , cex     = 0.8 * cex.ptsz
                            )#end legend
                     #---------------------------------------------------------------------#




                     #----- Plotting window and grid. -------------------------------------#
                     par(mar=c(5.1,4.4,4.1,2.1))
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                     axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                     axis(side=2,las=1)
                     title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.main)
                     if (plotgrid){
                        abline(v=when4,h=axTicks(2),col=grid.colour,lwd=0.75,lty="dotted")
                     }#end if
                     box()
                     #----- Plot the confidence interval. ---------------------------------#
                     for (r in sequence(nrate+1)){
                        do.call(what="epolygon",args=mci[[r]])
                     }#end for
                     #----- Plot the expected value. --------------------------------------#
                     for (r in sequence(nrate+1)){
                        do.call(what="points"  ,args=mdot[[r]])
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
                     #    Find the times for x scale.                                      #
                     #---------------------------------------------------------------------#
                     xlimit   = when.limit
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
                        legend ( x       = "bottom"
                               , inset   = 0.0
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , border  = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , col     = c(col.sta[1],rep("transparent",times=nrate-1)
                                            ,col.ed2[,1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(90,rep(0,times=nrate-1),angle)
                               , density = c(40,rep(0,times=nrate-1),dens)
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = cex.ptsz
                               )#end legend
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #    Loop over all DBH classes.                                    #
                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           left    = lo.box$left  [d ]
                           bottom  = lo.box$bottom[d ]
                           mar.now = lo.box$mar   [d,]





                           #----- Load the modelled rates. --------------------------------#
                           mult   = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                           mdot   = list()
                           mci    = list()
                           xlimit = when.limit
                           ylimit = NULL
                           ny     = n.census
                           for (r in sequence(nrate+1)){
                              if (r == nrate+1){
                                 mod.now      = sta[[sta.rate]][[indiv[i]]]$global
                                 cimult       = mult
                                 col.exp      = col.sta[1]
                                 col.ci       = col.sta[2]
                                 angle.ci     = 90
                                 dens.ci      = 40
                              }else{
                                 mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                                 cimult       = ifelse(ed22.ci,mult,NA)
                                 col.exp      = col.ed2[r,1]
                                 col.ci       = col.ed2[r,2]
                                 angle.ci     = angle[r]
                                 dens.ci      = dens [r]
                              }#end if

                              #----- Store polygons. --------------------------------------#
                              mdot[[r]] = list( x      = whenmid4
                                              , y      = mult*mod.now[1,d,-1]
                                              , col    = col.exp
                                              , pch    = 16
                                              , type   = "o"
                                              )#end list
                              mci [[r]] = list( x = c(whenmid4,rev(whenmid4))
                                              , y = cimult * c( mod.now[2,d,-1]
                                                              , rev(mod.now[3,d,-1])
                                                              )#end c
                                              , col   = col.ci
                                              , angle = angle.ci
                                              , dens  = dens.ci
                                              , lty   = "solid"
                                              , lwd   = 1.0
                                              )#end list
                              #------------------------------------------------------------#


                              #------------------------------------------------------------#
                              #     Update y range.                                        #
                              #------------------------------------------------------------#
                              if (global.ylim && r != (nrate + 1)){
                                 ylt.now = ylt[[ed2.rate[r]]][[indiv[i]]]
                                 ylt.min = mult * ylt.now$global$min[d]
                                 ylt.max = mult * ylt.now$global$max[d]
                                 ylimit  = range(c(ylimit,ylt.min,ylt.max))
                              }else{
                                 ylimit  = range( c(ylimit,mdot[[r]]$y,mci[[r]]$y)
                                                , na.rm=TRUE)
                              }#end if
                              #------------------------------------------------------------#
                           }#end for (r in sequence(nrate+1))
                           #------------------------------------------------------------------------#



                           #----- Set up the title and axes labels. -----------------------#
                           lesub = paste("DBH class:",dbh.names[d],sep="")
                           #---------------------------------------------------------------#


                           #----- Plot the box plot. --------------------------------------#
                           par(mar=c(4.1,3.1,2.1,2.1))
                           #----- Plotting window and grid. -------------------------------#
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                           axis(side=2,las=1)
                           title(main=lesub,xlab="",ylab="")
                           if (plotgrid){
                              abline(v=when4,h=axTicks(2),col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if
                           box()
                           #----- Plot the confidence interval. ---------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="epolygon",args=mci[[r]])
                           }#end for
                           #----- Plot the expected value. --------------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="points"  ,args=mdot[[r]])
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (d in 1:n.dbh)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
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
                              , cex.main  = 0.9*cex.ptsz
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




                  #========================================================================#
                  #========================================================================#
                  #    PFT-specific plots.                                                 #
                  #------------------------------------------------------------------------#
                  if (pfttoo){

                     #=====================================================================#
                     #=====================================================================#
                     #     PLOT-LEVEL rates.                                               #
                     #---------------------------------------------------------------------#
                     cat("  - PFT-dependent plot level: ",desc.indiv[i],"...","\n")
                     ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
                     sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(pftplot,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.pft = pretty.box(n=npfts+1,horizontal=TRUE)
                     #---------------------------------------------------------------------#





                     #----- Load the modelled rates. --------------------------------------#
                     mult   = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                     mdot   = list()
                     mci    = list()
                     xlimit = when.limit
                     ylimit = NULL
                     ny     = n.census
                     for (p in sequence(npfts+1)){
                        mdot[[p]] = list()
                        mci [[p]] = list()
                        for (r in sequence(nrate+1)){
                           if (r == nrate+1){
                              mod.now      = sta[[sta.rate]][[indiv[i]]]
                              cimult       = mult
                              col.exp      = col.sta[1]
                              col.ci       = col.sta[2]
                              angle.ci     = 90
                              dens.ci      = 40
                           }else{
                              mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]
                              cimult       = ifelse(ed22.ci,mult,NA)
                              col.exp      = col.ed2[r,1]
                              col.ci       = col.ed2[r,2]
                              angle.ci     = angle[r]
                              dens.ci      = dens [r]
                           }#end if

                           #----- Select the current PFT. ---------------------------------#
                           if (p == npfts+1){
                             expected.now = mult   * mod.now$global[1,-1]
                             q025.now     = cimult * mod.now$global[2,-1]
                             q975.now     = cimult * mod.now$global[3,-1]
                           }else{
                             expected.now = mult   * mod.now$expected[p,-1]
                             q025.now     = cimult * mod.now$q025    [p,-1]
                             q975.now     = cimult * mod.now$q975    [p,-1]
                           }#end if
                           #---------------------------------------------------------------#



                           #----- Store polygons. -----------------------------------------#
                           mdot[[p]][[r]] = list( x      = whenmid4
                                                , y      = expected.now
                                                , col    = col.exp
                                                , pch    = 16
                                                , type   = "o"
                                                )#end list
                           mci [[p]][[r]] = list( x = c(whenmid4,rev(whenmid4))
                                                , y = c(q025.now,rev(q975.now))
                                                , col   = col.ci
                                                , angle = angle.ci
                                                , dens  = dens.ci
                                                , lty   = "solid"
                                                , lwd   = 1.0
                                                )#end list
                           #---------------------------------------------------------------#


                           #---------------------------------------------------------------#
                           #     Update y range.                                           #
                           #---------------------------------------------------------------#
                           ylimit  = range( c(ylimit,mdot[[p]][[r]]$y,mci[[p]][[r]]$y)
                                          , na.rm=TRUE
                                          )#end range
                           #---------------------------------------------------------------#
                        }#end for (r in sequence(nrate+1))
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts+1))
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Loop over all formats, and make the plots.                      #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Open the file or the plot window. --------------------------#
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
                           pdf(file=fichier,onefile=FALSE,width=size$width
                              ,height=size$height,pointsize=ptsz,paper=size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Split the window.                                            #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par.orig = par(no.readonly = TRUE)
                        par(oma = c(0.2,3,4,0))
                        layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                              ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#


                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , border  = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , col     = c(col.sta[1],rep("transparent",times=nrate-1)
                                            ,col.ed2[,1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(90,rep(0,times=nrate-1),angle)
                               , density = c(40,rep(0,times=nrate-1),dens)
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 0.8 * cex.ptsz
                               )#end legend
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Loop over all PFTS to be plotted.                            #
                        #------------------------------------------------------------------#
                        for (p in sequence(npfts+1)){
                           if (p == npfts + 1){
                              pft.label = "All PFTs"
                           }else{
                              pft.label = pft$name[mypfts[p]]
                           }#end if
                           #---------------------------------------------------------------#


                           #---------------------------------------------------------------#
                           #     Update y range.                                           #
                           #---------------------------------------------------------------#
                           if (global.ylim){
                              ylimit = NULL
                              for (r in sequence(nrate)){
                                 ylt.now = ylt[[ed2.rate [r]]][[indiv[i]]]
                                 if (p == npfts+1){
                                    ylt.min = mult * ylt.now$global$min
                                    ylt.max = mult * ylt.now$global$max
                                 }else{
                                    ylt.min = mult * ylt.now$taxon$min[p]
                                    ylt.max = mult * ylt.now$taxon$max[p]
                                 }#end if
                                 #---------------------------------------------------------#
                                 ylimit  = range(c(ylimit,ylt.min,ylt.max))
                              }#end for
                              #------------------------------------------------------------#






                              #------------------------------------------------------------#
                              #     Define margins.                                        #
                              #------------------------------------------------------------#
                              left    = TRUE
                              bottom  = TRUE
                              mar.now = lo.pft$mar0
                              #------------------------------------------------------------#
                           }else{
                              #------------------------------------------------------------#
                              #     Define margins.                                        #
                              #------------------------------------------------------------#
                              left    = lo.pft$left  [p ]
                              bottom  = lo.pft$bottom[p ]
                              mar.now = lo.pft$mar   [p,]
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#




                           #----- Plotting window and grid. -------------------------------#
                           par(mar=mar.now)
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           if (plotgrid){
                              abline(v=when4,h=axTicks(2),col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if
                           if (bottom) axis.rt(side=1,at=when4,labels=when.label
                                              ,las=5,off=.05)
                           if (left)   axis(side=2,las=1)
                           title(main=pft.label)
                           box()
                           #----- Plot the confidence interval. ---------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="epolygon",args=mci[[p]][[r]])
                           }#end for
                           #----- Plot the expected value. --------------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="points"  ,args=mdot[[p]][[r]])
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (p in sequence(npfts)
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
                        #------------------------------------------------------------------#


                        #------------------------------------------------------------------#
                        #      Plot global annotation.                                     #
                        #------------------------------------------------------------------#
                        gtitle( main      = letitre
                              , xlab      = lex
                              , ylab      = ley
                              , off.xlab  = 1/6
                              , line.xlab = 4.1
                              , line.ylab = 2.6
                              , cex.main  = 0.9*cex.ptsz
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
                     }#end for (o in 1:nout)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo)
                  #========================================================================#
                  #========================================================================#



                  #========================================================================#
                  #========================================================================#
                  #     DBH-LEVEL rates by PFT.                                            #
                  #------------------------------------------------------------------------#
                  if (sizetoo && pfttoo){
                     cat("  - PFT-dependent, DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(pftsize,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Loop over all PFTs.                                              #
                     #---------------------------------------------------------------------#
                     for (p in sequence(npfts)){
                        pft.label = pft$name[mypfts[p]]
                        pft.key   = paste("pft",sprintf("%2.2i",mypfts[p]),sep="")



                        #------------------------------------------------------------------#
                        #     Loop over all formats, and make the plots.                   #
                        #------------------------------------------------------------------#
                        for (o in 1:nout){
                           #----- Open the file or the plot window. -----------------------#
                           fichier = file.path(outindiv
                                              ,paste(iata,"-theme-",theme.now,"-",indiv[i],
                                                    "-",pft.key,"-",iphen.key[ph]
                                                    ,".",outform[o] ,sep=""))
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
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Split the window into several smaller windows.            #
                           #---------------------------------------------------------------#
                           par(par.user)
                           par.orig = par(no.readonly = TRUE)
                           par(oma = c(0.2,3,4,0))
                           layout(mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
                                 ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                                 )#end layout
                           #---------------------------------------------------------------#


                           #----- Plot legend. --------------------------------------------#
                           if (ed22.ci){
                              par(mar=c(0.1,0.1,0.1,0.1))
                              plot.new()
                              plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                              legend ( x       = "bottom"
                                     , inset   = 0.0
                                     , legend  = c("Observed",rep("",times=nrate-1)
                                                  ,desc.rate)
                                     , fill    = c(col.sta[2]
                                                  ,rep("transparent",times=nrate-1)
                                                  ,col.ed2[,2])
                                     , border  = c(col.sta[2]
                                                  ,rep("transparent",times=nrate-1)
                                                  ,col.ed2[,2])
                                     , col     = c(col.sta[1]
                                                  ,rep("transparent",times=nrate-1)
                                                  ,col.ed2[,1])
                                     , lwd     = 2.0
                                     , pt.cex  = 1.0
                                     , angle   = c(90,rep(0,times=nrate-1),angle)
                                     , density = c(40,rep(0,times=nrate-1),dens)
                                     , ncol    = 2
                                     , title   = "(Shaded - 95% C.I.)"
                                     , cex     = 0.8 * cex.ptsz
                                     )#end legend
                           }else{
                              par(mar=c(0.1,0.1,0.1,0.1))
                              plot.new()
                              plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                              legend ( x       = "bottom"
                                     , inset   = 0.0
                                     , legend  = c("Observed",rep("",times=nrate-1)
                                                  ,desc.rate)
                                     , fill    = c(col.sta[2]
                                                  ,rep("transparent",times=nrate+1))
                                     , border  = c(col.sta[2]
                                                  ,rep("transparent",times=nrate+1))
                                     , col     = c(col.sta[1]
                                                  ,rep("transparent",times=nrate-1))
                                     , lwd     = 2.0
                                     , pt.cex  = 1.0
                                     , angle   = c(90,rep(0,times=nrate+1))
                                     , density = c(40,rep(0,times=nrate+1))
                                     , ncol    = 2
                                     , title   = "(Shaded - 95% C.I.)"
                                     , cex     = 0.8 * cex.ptsz
                                     )#end legend
                           }#end if
                           #---------------------------------------------------------------#




                           #---------------------------------------------------------------#
                           #    Loop over all DBH classes.                                 #
                           #---------------------------------------------------------------#
                           for (d in 1:n.dbh){

                              #----- Load the modelled rates. -----------------------------#
                              mult = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                              mdot = list()
                              mci  = list()
                              xlimit = when.limit
                              ylimit = NULL
                              ny     = n.census
                              for (r in sequence(nrate+1)){
                                 if (r == nrate+1){
                                    mod.now      = sta[[sta.rate]][[indiv[i]]]
                                    cimult       = mult
                                    col.exp      = col.sta[1]
                                    col.ci       = col.sta[2]
                                    angle.ci     = 90
                                    dens.ci      = 40
                                 }else{
                                    mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]
                                    cimult       = ifelse(ed22.ci,mult,NA)
                                    col.exp      = col.ed2[r,1]
                                    col.ci       = col.ed2[r,2]
                                    angle.ci     = angle[r]
                                    dens.ci      = dens [r]
                                 }#end if

                                 #----- Select the current PFT. ---------------------------#
                                 expected.now = mult   * mod.now$expected[p,d,-1]
                                 q025.now     = cimult * mod.now$q025    [p,d,-1]
                                 q975.now     = cimult * mod.now$q975    [p,d,-1]
                                 #---------------------------------------------------------#



                                 #----- Store polygons. -----------------------------------#
                                 mdot[[r]] = list( x      = whenmid4
                                                 , y      = expected.now
                                                 , col    = col.exp
                                                 , pch    = 16
                                                 , type   = "o"
                                                 )#end list
                                 mci [[r]] = list( x     = c(whenmid4,rev(whenmid4))
                                                 , y     = c(q025.now,rev(q975.now))
                                                 , col   = col.ci
                                                 , angle = angle.ci
                                                 , dens  = dens.ci
                                                 , lty   = "solid"
                                                 , lwd   = 1.0
                                                 )#end list
                                 #---------------------------------------------------------#


                                 #---------------------------------------------------------------------#
                                 #     Update y range.                                                 #
                                 #---------------------------------------------------------------------#
                                 if (global.ylim && r != (nrate + 1)){
                                    ylt.now = ylt[[ed2.rate[r]]][[indiv[i]]]
                                    ylt.min = mult * ylt.now$taxon$min[p,d]
                                    ylt.max = mult * ylt.now$taxon$max[p,d]
                                    ylimit  = range(c(ylimit,ylt.min,ylt.max))
                                 }else{
                                    ylimit  = range( c(ylimit,mdot[[r]]$y,mci[[r]]$y)
                                                   , na.rm = TRUE
                                                   )#end range
                                 }#end if
                                 #---------------------------------------------------------------------#
                              }#end for (r in sequence(nrate+1))
                              #------------------------------------------------------------#


                              #----- Set up the title and axes labels. --------------------#
                              lesub = paste("DBH class:",dbh.names[d],sep="")
                              #------------------------------------------------------------#


                              #----- Plotting window and grid. ----------------------------#
                              par(mar=lo.box$mar0)
                              plot.new()
                              plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                              if (plotgrid){
                                 abline(v=when4,h=axTicks(2),col=grid.colour
                                       ,lwd=0.75,lty="dotted")
                              }#end if
                              axis.rt(side=1,at=when4,labels=when.label,las=5
                                     ,off=.05)
                              axis(side=2,las=1)
                              box()
                              title(main=lesub,xlab="",ylab="")
                              #----- Plot the confidence interval. ------------------------#
                              for (r in sequence(nrate+1)){
                                 do.call(what="epolygon",args=mci[[r]])
                              }#end for
                              #----- Plot the expected value. -----------------------------#
                              for (r in sequence(nrate+1)){
                                 do.call(what="points"  ,args=mdot[[r]])
                              }#end for
                              #------------------------------------------------------------#
                           }#end for (d in 1:n.dbh)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Make the title and axis labels.                           #
                           #---------------------------------------------------------------#
                           letitre = paste(theme.desc," - ",longname," - ",pft.label
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                           ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                           lex     = "" # desc.unit(desc="Census year",unit=NULL)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Plot title.                                               #
                           #---------------------------------------------------------------#
                           gtitle( main      = letitre
                                 , xlab      = lex
                                 , ylab      = ley
                                 , off.xlab  = 1/6
                                 , line.xlab = 4.1
                                 , line.ylab = 2.6
                                 , cex.main  = 0.9*cex.ptsz
                                 )#end gtitle
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for (o in 1:nout)
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo && sizetoo)
                  #------------------------------------------------------------------------#
               }#end for (i in 1:nindiv)
               #---------------------------------------------------------------------------#
            }#end for
            #==============================================================================#
            #==============================================================================#





            #==============================================================================#
            #==============================================================================#
            #  8. Plot the theme time series of the rates (by DBH class if applicable).    #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratetheme){
               this.plot  = pratetheme[[n]]
               sizetoo    = this.plot$sizetoo
               pfttoo     = this.plot$pfttoo
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
               col.ed2    = matrix(col.ed2,nrow=nrate)

               cat(" + Plotting flat time series of ",theme.desc,"...","\n")
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
                  mult   = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                  mexp   = list()
                  mdot   = list()
                  mci    = list()
                  xlimit = when.limit
                  ylimit = NULL
                  ny     = n.census
                  for (r in sequence(nrate+1)){
                     if (r == nrate+1){
                        mod.now      = sta[[sta.rate]][[indiv[i]]]$global
                        cimult       = mult
                        col.exp      = col.sta[1]
                        col.ci       = col.sta[2]
                        angle.ci     = 90
                        dens.ci      = 40
                     }else{
                        mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                        cimult       = ifelse(ed22.ci,mult,NA)
                        col.exp      = col.ed2[r,1]
                        col.ci       = col.ed2[r,2]
                        angle.ci     = angle[r]
                        dens.ci      = dens [r]
                     }#end if

                     #----- Store polygons. -----------------------------------------------#
                     mexp[[r]] = list( x0     = when4[-ny]
                                     , x1     = when4[ -1]
                                     , y0     = mult*mod.now[1,-1]
                                     , y1     = mult*mod.now[1,-1]
                                     , col    = col.exp
                                     , lwd    = 2.0
                                     )#end list
                     mdot[[r]] = list( x      = whenmid4
                                     , y      = mult*mod.now[1,-1]
                                     , col    = col.exp
                                     , pch    = 16
                                     , type   = "p"
                                     )#end list
                     mci [[r]] = list( x = c( rbind( when4[-ny], when4[ -1]
                                                   , when4[ -1], when4[-ny]
                                                   , NA) )
                                     , y = cimult * c( rbind( mod.now[2,-1],mod.now[2,-1]
                                                            , mod.now[3,-1],mod.now[3,-1]
                                                            , NA )
                                                     )#end c
                                     , col   = col.ci
                                     , angle = angle.ci
                                     , dens  = dens.ci
                                     , lty   = "solid"
                                     , lwd   = 1.0
                                     )#end list
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Update y range.                                                 #
                     #---------------------------------------------------------------------#
                     if (global.ylim && r != (nrate + 1)){
                        ylt.min = mult * ylt[[ed2.rate[r]]][[indiv[i]]]$global$min
                        ylt.max = mult * ylt[[ed2.rate[r]]][[indiv[i]]]$global$max
                        ylimit  = range(c(ylimit,ylt.min,ylt.max))
                     }else{
                        ylimit  = range(c(ylimit,mdot[[r]]$y,mci[[r]]$y),na.rm=TRUE)
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for (r in sequence(nrate+1))
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Loop over all formats, and make the plots.                         #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Open the file or the plot window. -----------------------------#
                     fichier = file.path(outindiv
                                        ,paste(iata,"-thflt-",theme.now,"-",indiv[i]
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
                     letitre = paste(theme.desc," - ",longname
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                     ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                     lex     = "" # desc.unit(desc="Census year",unit=NULL)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Go on and plot stuff.                                           #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     layout(mat=rbind(2,1),heights=c(5,1))
                     #---------------------------------------------------------------------#


                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.1,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                            , fill    = c(col.sta[2],rep("transparent",times=nrate-1)
                                         ,col.ed2[,2])
                            , border  = c(col.sta[2],rep("transparent",times=nrate-1)
                                         ,col.ed2[,2])
                            , col     = c(col.sta[1],rep("transparent",times=nrate-1)
                                         ,col.ed2[,1])
                            , lwd     = 2.0
                            , pt.cex  = 1.0
                            , angle   = c(90,rep(0,times=nrate-1),angle)
                            , density = c(40,rep(0,times=nrate-1),dens)
                            , ncol    = 2
                            , title   = "(Shaded - 95% C.I.)"
                            , cex     = 0.8 * cex.ptsz
                            )#end legend
                     #---------------------------------------------------------------------#




                     #----- Plotting window and grid. -------------------------------------#
                     par(mar=c(5.1,4.4,4.1,2.1))
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                     axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                     axis(side=2,las=1)
                     title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.main)
                     if (plotgrid){
                        abline(v=when4,h=axTicks(2),col=grid.colour,lwd=0.75,lty="dotted")
                     }#end if
                     box()
                     #----- Plot the confidence interval. ---------------------------------#
                     for (r in sequence(nrate+1)){
                        do.call(what="epolygon",args=mci[[r]])
                     }#end for
                     #----- Plot the expected value. --------------------------------------#
                     for (r in sequence(nrate+1)){
                        do.call(what="segments",args=mexp[[r]])
                        do.call(what="points"  ,args=mdot[[r]])
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
                     xlimit   = when.limit
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
                                           ,paste(iata,"-thflt-",theme.now,"-",indiv[i],"-"
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
                        legend ( x       = "bottom"
                               , inset   = 0.0
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , border  = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , col     = c(col.sta[1],rep("transparent",times=nrate-1)
                                            ,col.ed2[,1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(90,rep(0,times=nrate-1),angle)
                               , density = c(40,rep(0,times=nrate-1),dens)
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = cex.ptsz
                               )#end legend
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #    Loop over all DBH classes.                                    #
                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){
                           #----- Load the modelled rates. --------------------------------#
                           mult   = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                           mexp   = list()
                           mdot   = list()
                           mci    = list()
                           xlimit = when.limit
                           ylimit = NULL
                           ny     = n.census
                           for (r in sequence(nrate+1)){
                              if (r == nrate+1){
                                 mod.now      = sta[[sta.rate]][[indiv[i]]]$global
                                 cimult       = mult
                                 col.exp      = col.sta[1]
                                 col.ci       = col.sta[2]
                                 angle.ci     = 90
                                 dens.ci      = 40
                              }else{
                                 mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                                 cimult       = ifelse(ed22.ci,mult,NA)
                                 col.exp      = col.ed2[r,1]
                                 col.ci       = col.ed2[r,2]
                                 angle.ci     = angle[r]
                                 dens.ci      = dens [r]
                              }#end if

                              #----- Store polygons. --------------------------------------#
                              mexp[[r]] = list( x0     = when4[-ny]
                                              , x1     = when4[ -1]
                                              , y0     = mult*mod.now[1,d,-1]
                                              , y1     = mult*mod.now[1,d,-1]
                                              , col    = col.exp
                                              , lwd    = 2.0
                                              )#end list
                              mdot[[r]] = list( x      = whenmid4
                                              , y      = mult*mod.now[1,d,-1]
                                              , col    = col.exp
                                              , pch    = 16
                                              , type   = "p"
                                              )#end list
                              mci [[r]] = list( x = c( rbind( when4[-ny], when4[ -1]
                                                            , when4[ -1], when4[-ny]
                                                            , NA) )
                                              , y = cimult * c( rbind( mod.now[2,d,-1]
                                                                     , mod.now[2,d,-1]
                                                                     , mod.now[3,d,-1]
                                                                     , mod.now[3,d,-1]
                                                                     , NA )
                                                              )#end c
                                              , col   = col.ci
                                              , angle = angle.ci
                                              , dens  = dens.ci
                                              , lty   = "solid"
                                              , lwd   = 1.0
                                              )#end list
                              #------------------------------------------------------------#


                              #------------------------------------------------------------#
                              #     Update y range.                                        #
                              #------------------------------------------------------------#
                              if (global.ylim && r != (nrate + 1)){
                                 ylt.now = ylt[[ed2.rate[r]]][[indiv[i]]]
                                 ylt.min = mult * ylt.now$global$min[d]
                                 ylt.max = mult * ylt.now$global$max[d]
                                 ylimit  = range(c(ylimit,ylt.min,ylt.max))
                              }else{
                                 ylimit  = range( c(ylimit,mdot[[r]]$y,mci[[r]]$y)
                                                , na.rm=TRUE)
                              }#end if
                              #------------------------------------------------------------#
                           }#end for (r in sequence(nrate+1))
                           #------------------------------------------------------------------------#



                           #----- Set up the title and axes labels. -----------------------#
                           lesub = paste("DBH class:",dbh.names[d],sep="")
                           #---------------------------------------------------------------#


                           #----- Plotting window and grid. -------------------------------#
                           par(mar=lo.box$mar0)
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                           axis(side=2,las=1)
                           title(main=lesub,xlab="",ylab="")
                           if (plotgrid){
                              abline(v=when4,h=axTicks(2),col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if
                           box()
                           #----- Plot the confidence interval. ---------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="epolygon",args=mci[[r]])
                           }#end for
                           #----- Plot the expected value. --------------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="segments",args=mexp[[r]])
                              do.call(what="points"  ,args=mdot[[r]])
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (d in 1:n.dbh)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
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
                              , cex.main  = 0.9*cex.ptsz
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





                  #========================================================================#
                  #========================================================================#
                  #    PFT-specific plots.                                                 #
                  #------------------------------------------------------------------------#
                  if (pfttoo){

                     #=====================================================================#
                     #=====================================================================#
                     #     PLOT-LEVEL rates.                                               #
                     #---------------------------------------------------------------------#
                     cat("  - PFT-dependent plot level: ",desc.indiv[i],"...","\n")
                     ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
                     sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(pftplot,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.pft = pretty.box(n=npfts+1,horizontal=TRUE)
                     #---------------------------------------------------------------------#





                     #----- Load the modelled rates. --------------------------------------#
                     mult   = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                     mexp   = list()
                     mdot   = list()
                     mci    = list()
                     xlimit = when.limit
                     ylimit = NULL
                     ny     = n.census
                     for (p in sequence(npfts+1)){
                        mexp[[p]] = list()
                        mdot[[p]] = list()
                        mci [[p]] = list()
                        for (r in sequence(nrate+1)){
                           if (r == nrate+1){
                              mod.now      = sta[[sta.rate]][[indiv[i]]]
                              cimult       = mult
                              col.exp      = col.sta[1]
                              col.ci       = col.sta[2]
                              angle.ci     = 90
                              dens.ci      = 40
                           }else{
                              mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]
                              cimult       = ifelse(ed22.ci,mult,NA)
                              col.exp      = col.ed2[r,1]
                              col.ci       = col.ed2[r,2]
                              angle.ci     = angle[r]
                              dens.ci      = dens [r]
                           }#end if

                           #----- Select the current PFT. ---------------------------------#
                           if (p == npfts+1){
                             expected.now = mult   * mod.now$global[1,-1]
                             q025.now     = cimult * mod.now$global[2,-1]
                             q975.now     = cimult * mod.now$global[3,-1]
                           }else{
                             expected.now = mult   * mod.now$expected[p,-1]
                             q025.now     = cimult * mod.now$q025    [p,-1]
                             q975.now     = cimult * mod.now$q975    [p,-1]
                           }#end if
                           #---------------------------------------------------------------#



                           #----- Store polygons. -----------------------------------------#
                           mexp[[p]][[r]] = list( x0     = when4[-ny]
                                                , x1     = when4[ -1]
                                                , y0     = expected.now
                                                , y1     = expected.now
                                                , col    = col.exp
                                                , lwd    = 2.0
                                                )#end list
                           mdot[[p]][[r]] = list( x      = whenmid4
                                                , y      = expected.now
                                                , col    = col.exp
                                                , pch    = 16
                                                , type   = "p"
                                                )#end list
                           mci [[p]][[r]] = list( x = c( rbind( when4[-ny], when4[ -1]
                                                              , when4[ -1], when4[-ny]
                                                              , NA) )
                                                , y = c( rbind( q025.now,q025.now
                                                              , q975.now,q975.now
                                                              , NA )
                                                       )#end c
                                                , col   = col.ci
                                                , angle = angle.ci
                                                , dens  = dens.ci
                                                , lty   = "solid"
                                                , lwd   = 1.0
                                                )#end list
                           #---------------------------------------------------------------#


                           #---------------------------------------------------------------#
                           #     Update y range.                                           #
                           #---------------------------------------------------------------#
                           ylimit  = range( c(ylimit,mdot[[p]][[r]]$y,mci[[p]][[r]]$y)
                                          , na.rm=TRUE
                                          )#end range
                           #---------------------------------------------------------------#
                        }#end for (r in sequence(nrate+1))
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts+1))
                     #---------------------------------------------------------------------#




                     #---------------------------------------------------------------------#
                     #     Loop over all formats, and make the plots.                      #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Open the file or the plot window. --------------------------#
                        fichier = file.path(outindiv
                                           ,paste(iata,"-thflt-",theme.now,"-",indiv[i]
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
                           pdf(file=fichier,onefile=FALSE,width=size$width
                              ,height=size$height,pointsize=ptsz,paper=size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Split the window.                                            #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par.orig = par(no.readonly = TRUE)
                        par(oma = c(0.2,3,4,0))
                        layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                              ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#


                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , border  = c(col.sta[2],rep("transparent",times=nrate-1)
                                            ,col.ed2[,2])
                               , col     = c(col.sta[1],rep("transparent",times=nrate-1)
                                            ,col.ed2[,1])
                               , lwd     = 2.0
                               , pt.cex  = 1.0
                               , angle   = c(90,rep(0,times=nrate-1),angle)
                               , density = c(40,rep(0,times=nrate-1),dens)
                               , ncol    = 2
                               , title   = "(Shaded - 95% C.I.)"
                               , cex     = 0.8 * cex.ptsz
                               )#end legend
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Loop over all PFTS to be plotted.                            #
                        #------------------------------------------------------------------#
                        for (p in sequence(npfts+1)){
                           if (p == npfts + 1){
                              pft.label = "All PFTs"
                           }else{
                              pft.label = pft$name[mypfts[p]]
                           }#end if


                           #---------------------------------------------------------------#
                           #     Update y range.                                           #
                           #---------------------------------------------------------------#
                           if (global.ylim){
                              ylimit = NULL
                              for (r in sequence(nrate)){
                                 ylt.now = ylt[[ed2.rate [r]]][[indiv[i]]]
                                 if (p == npfts+1){
                                    ylt.min = mult * ylt.now$global$min
                                    ylt.max = mult * ylt.now$global$max
                                 }else{
                                    ylt.min = mult * ylt.now$taxon$min[p]
                                    ylt.max = mult * ylt.now$taxon$max[p]
                                 }#end if
                                 #---------------------------------------------------------#
                                 ylimit  = range(c(ylimit,ylt.min,ylt.max))
                              }#end for
                              #------------------------------------------------------------#
                              #------------------------------------------------------------#
                              #     Define margins.                                        #
                              #------------------------------------------------------------#
                              left    = TRUE
                              bottom  = TRUE
                              mar.now = lo.pft$mar0
                              #------------------------------------------------------------#
                           }else{


                              #------------------------------------------------------------#
                              #     Define margins.                                        #
                              #------------------------------------------------------------#
                              left    = lo.pft$left  [p ]
                              bottom  = lo.pft$bottom[p ]
                              mar.now = lo.pft$mar   [p,]
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#




                           #----- Plotting window and grid. -------------------------------#
                           par(mar=mar.now)
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           if (plotgrid){
                              abline(v=when4,h=axTicks(2),col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if
                           if (bottom) axis.rt(side=1,at=when4,labels=when.label
                                              ,las=5,off=.05)
                           if (left)   axis(side=2,las=1)
                           title(main=pft.label)
                           box()
                           #----- Plot the confidence interval. ---------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="epolygon",args=mci[[p]][[r]])
                           }#end for
                           #----- Plot the expected value. --------------------------------#
                           for (r in sequence(nrate+1)){
                              do.call(what="segments",args=mexp[[p]][[r]])
                              do.call(what="points"  ,args=mdot[[p]][[r]])
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (p in sequence(npfts)
                        #------------------------------------------------------------------#






                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
                        #------------------------------------------------------------------#


                        #------------------------------------------------------------------#
                        #      Plot global annotation.                                     #
                        #------------------------------------------------------------------#
                        gtitle( main      = letitre
                              , xlab      = lex
                              , ylab      = ley
                              , off.xlab  = 1/6
                              , line.xlab = 4.1
                              , line.ylab = 2.6
                              , cex.main  = 0.9*cex.ptsz
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
                     }#end for (o in 1:nout)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo)
                  #========================================================================#
                  #========================================================================#




                  #========================================================================#
                  #========================================================================#
                  #     DBH-LEVEL rates by PFT.                                            #
                  #------------------------------------------------------------------------#
                  if (sizetoo && pfttoo){
                     cat("  - PFT-dependent, DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(pftsize,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Loop over all PFTs.                                              #
                     #---------------------------------------------------------------------#
                     for (p in sequence(npfts)){
                        pft.label = pft$name[mypfts[p]]
                        pft.key   = paste("pft",sprintf("%2.2i",mypfts[p]),sep="")



                        #------------------------------------------------------------------#
                        #     Loop over all formats, and make the plots.                   #
                        #------------------------------------------------------------------#
                        for (o in 1:nout){
                           #----- Open the file or the plot window. -----------------------#
                           fichier = file.path(outindiv
                                              ,paste(iata,"-thflt-",theme.now,"-",indiv[i],
                                                    "-",pft.key,"-",iphen.key[ph]
                                                    ,".",outform[o] ,sep=""))
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
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Split the window into several smaller windows.            #
                           #---------------------------------------------------------------#
                           par(par.user)
                           par.orig = par(no.readonly = TRUE)
                           par(oma = c(0.2,3,4,0))
                           layout(mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
                                 ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                                 )#end layout
                           #---------------------------------------------------------------#


                           #----- Plot legend. --------------------------------------------#
                           if (ed22.ci){
                              par(mar=c(0.1,0.1,0.1,0.1))
                              plot.new()
                              plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                              legend ( x       = "bottom"
                                     , inset   = 0.0
                                     , legend  = c("Observed",rep("",times=nrate-1)
                                                  ,desc.rate)
                                     , fill    = c(col.sta[2]
                                                  ,rep("transparent",times=nrate-1)
                                                  ,col.ed2[,2])
                                     , border  = c(col.sta[2]
                                                  ,rep("transparent",times=nrate-1)
                                                  ,col.ed2[,2])
                                     , col     = c(col.sta[1]
                                                  ,rep("transparent",times=nrate-1)
                                                  ,col.ed2[,1])
                                     , lwd     = 2.0
                                     , pt.cex  = 1.0
                                     , angle   = c(90,rep(0,times=nrate-1),angle)
                                     , density = c(40,rep(0,times=nrate-1),dens)
                                     , ncol    = 2
                                     , title   = "(Shaded - 95% C.I.)"
                                     , cex     = 0.8 * cex.ptsz
                                     )#end legend
                           }else{
                              par(mar=c(0.1,0.1,0.1,0.1))
                              plot.new()
                              plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                              legend ( x       = "bottom"
                                     , inset   = 0.0
                                     , legend  = c("Observed",rep("",times=nrate-1)
                                                  ,desc.rate)
                                     , fill    = c(col.sta[2]
                                                  ,rep("transparent",times=nrate+1))
                                     , border  = c(col.sta[2]
                                                  ,rep("transparent",times=nrate+1))
                                     , col     = c(col.sta[1]
                                                  ,rep("transparent",times=nrate-1))
                                     , lwd     = 2.0
                                     , pt.cex  = 1.0
                                     , angle   = c(90,rep(0,times=nrate+1))
                                     , density = c(40,rep(0,times=nrate+1))
                                     , ncol    = 2
                                     , title   = "(Shaded - 95% C.I.)"
                                     , cex     = 0.8 * cex.ptsz
                                     )#end legend
                           }#end if
                           #---------------------------------------------------------------#




                           #---------------------------------------------------------------#
                           #    Loop over all DBH classes.                                 #
                           #---------------------------------------------------------------#
                           for (d in 1:n.dbh){
                              #----- Load the modelled rates. -----------------------------#
                              mult = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                              mexp = list()
                              mdot = list()
                              mci  = list()
                              xlimit = when.limit
                              ylimit = NULL
                              ny     = n.census
                              for (r in sequence(nrate+1)){
                                 if (r == nrate+1){
                                    mod.now      = sta[[sta.rate]][[indiv[i]]]
                                    cimult       = mult
                                    col.exp      = col.sta[1]
                                    col.ci       = col.sta[2]
                                    angle.ci     = 90
                                    dens.ci      = 40
                                 }else{
                                    mod.now      = ed2[[ed2.rate[r]]][[indiv[i]]]
                                    cimult       = ifelse(ed22.ci,mult,NA)
                                    col.exp      = col.ed2[r,1]
                                    col.ci       = col.ed2[r,2]
                                    angle.ci     = angle[r]
                                    dens.ci      = dens [r]
                                 }#end if

                                 #----- Select the current PFT. ---------------------------#
                                 expected.now = mult   * mod.now$expected[p,d,-1]
                                 q025.now     = cimult * mod.now$q025    [p,d,-1]
                                 q975.now     = cimult * mod.now$q975    [p,d,-1]
                                 #---------------------------------------------------------#



                                 #----- Store polygons. -----------------------------------#
                                 mexp[[r]] = list( x0     = when4[-ny]
                                                 , x1     = when4[ -1]
                                                 , y0     = expected.now
                                                 , y1     = expected.now
                                                 , col    = col.exp
                                                 , lwd    = 2.0
                                                 )#end list
                                 mdot[[r]] = list( x      = whenmid4
                                                 , y      = expected.now
                                                 , col    = col.exp
                                                 , pch    = 16
                                                 , type   = "p"
                                                 )#end list
                                 mci [[r]] = list( x = c( rbind( when4[-ny], when4[ -1]
                                                               , when4[ -1], when4[-ny]
                                                               , NA) )
                                                 , y = c( rbind( q025.now,q025.now
                                                               , q975.now,q975.now
                                                               , NA )
                                                        )#end c
                                                 , col   = col.ci
                                                 , angle = angle.ci
                                                 , dens  = dens.ci
                                                 , lty   = "solid"
                                                 , lwd   = 1.0
                                                 )#end list
                                 #---------------------------------------------------------#


                                 #---------------------------------------------------------------------#
                                 #     Update y range.                                                 #
                                 #---------------------------------------------------------------------#
                                 if (global.ylim && r != (nrate + 1)){
                                    ylt.now = ylt[[ed2.rate[r]]][[indiv[i]]]
                                    ylt.min = mult * ylt.now$taxon$min[p,d]
                                    ylt.max = mult * ylt.now$taxon$max[p,d]
                                    ylimit  = range(c(ylimit,ylt.min,ylt.max))
                                 }else{
                                    ylimit  = range( c(ylimit,mdot[[r]]$y,mci[[r]]$y)
                                                   , na.rm = TRUE
                                                   )#end range
                                 }#end if
                                 #---------------------------------------------------------------------#
                              }#end for (r in sequence(nrate+1))
                              #------------------------------------------------------------#


                              #----- Set up the title and axes labels. --------------------#
                              lesub = paste("DBH class:",dbh.names[d],sep="")
                              #------------------------------------------------------------#


                              #----- Plot the box plot. -----------------------------------#
                              par(mar=lo.box$mar0)
                              #----- Plotting window and grid. ----------------------------#
                              plot.new()
                              plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                              if (plotgrid){
                                 abline(v=when4,h=axTicks(2),col=grid.colour
                                       ,lwd=0.75,lty="dotted")
                              }#end if
                              axis.rt(side=1,at=when4,labels=when.label
                                     ,las=5,off=.05)
                              axis(side=2,las=1)
                              box()
                              title(main=lesub,xlab="",ylab="")
                              #----- Plot the confidence interval. ------------------------#
                              for (r in sequence(nrate+1)){
                                 do.call(what="epolygon",args=mci[[r]])
                              }#end for
                              #----- Plot the expected value. -----------------------------#
                              for (r in sequence(nrate+1)){
                                 do.call(what="segments",args=mexp[[r]])
                                 do.call(what="points"  ,args=mdot[[r]])
                              }#end for
                              #------------------------------------------------------------#
                           }#end for (d in 1:n.dbh)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Make the title and axis labels.                           #
                           #---------------------------------------------------------------#
                           letitre = paste(theme.desc," - ",longname," - ",pft.label
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                           ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                           lex     = "" # desc.unit(desc="Census year",unit=NULL)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Plot title.                                               #
                           #---------------------------------------------------------------#
                           gtitle( main      = letitre
                                 , xlab      = lex
                                 , ylab      = ley
                                 , off.xlab  = 1/6
                                 , line.xlab = 4.1
                                 , line.ylab = 2.6
                                 , cex.main  = 0.9*cex.ptsz
                                 )#end gtitle
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for (o in 1:nout)
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo && sizetoo)
                  #------------------------------------------------------------------------#
               }#end for (i in 1:nindiv)
               #---------------------------------------------------------------------------#
            }#end for
            #==============================================================================#
            #==============================================================================#





            #==============================================================================#
            #==============================================================================#
            #  9. Plot the time series of the rates (by DBH class if applicable).          #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratets){
               this.plot  = pratets[[n]]
               pfttoo     = this.plot$pfttoo
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
                  mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                  sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                  sta.expected = mult * sta.mod[1,-1]
                  sta.q025     = mult * sta.mod[2,-1]
                  sta.q975     = mult * sta.mod[3,-1]
                  ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
                  ed2.expected = mult * ed2.mod[1,-1]
                  ed2.q025     = mult * ed2.mod[2,-1]
                  ed2.q975     = mult * ed2.mod[3,-1]
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #    Find the DBH for x scale.                                           #
                  #------------------------------------------------------------------------#
                  x.years  = whenmid4
                  xlimit   = when.limit
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
                     yuse   = c(sta.q025,sta.q975,sta.expected
                               ,ed2.q025,ed2.q975,ed2.expected)
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
                     yuse   = c(sta.q025,sta.q975,sta.expected,ed2.expected)
                     ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Over-write with global rate.                                       #
                  #------------------------------------------------------------------------#
                  if (global.ylim){
                     ylt.mod = ylt[[ed2.rate]][[indiv[i]]]$global
                     ylt.min = mult * ylt.mod$min
                     ylt.max = mult * ylt.mod$max
                     ylimit  = pretty.xylim(u=c(ylt.min,ylt.max),fracexp=0.0,is.log=ylog)
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
                     letitre = paste(desc.rate," - ",longname
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                     ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                     lex     = "" # desc.unit(desc="Census year",unit=NULL)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Go on and plot stuff.                                           #
                     #---------------------------------------------------------------------#
                     par(mar=c(5.1,4.4,4.1,2.1))
                     #----- Plotting window and grid. -------------------------------------#
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                     title(main=letitre,xlab=lex,ylab=ley,log=plog,cex.main=0.7*cex.main)
                     if (plotgrid){
                        abline(h=axTicks(2),v=when4,col=grid.colour,lwd=0.75,lty="dotted")
                     }#end if
                     axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                     axis(side=2,las=1)
                     #----- Plot the taxon rate with confidence interval. -----------------#
                     epolygon(x=plot.poly$x,y=plot.poly$y,col=plot.poly$col,angle=c(-45,45)
                             ,density=40,lty="solid",lwd=1.0)
                     lines (x=whenmid4,y=sta.expected,lwd=2.0,col=col.sta[1])
                     lines (x=whenmid4,y=ed2.expected,lwd=2.0,col=col.ed2[1])
                     points(x=whenmid4,y=sta.expected,pch=16 ,col=col.sta[1])
                     points(x=whenmid4,y=ed2.expected,pch=16 ,col=col.ed2[1])


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
                     mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                     sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                     sta.expected = mult * sta.mod[1,,-1]
                     sta.q025     = mult * sta.mod[2,,-1]
                     sta.q975     = mult * sta.mod[3,,-1]
                     ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
                     ed2.expected = mult * ed2.mod[1,,-1]
                     ed2.q025     = mult * ed2.mod[2,,-1]
                     ed2.q975     = mult * ed2.mod[3,,-1]
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Find the DBH for x scale.                                        #
                     #---------------------------------------------------------------------#
                     x.years  = whenmid4
                     xlimit   = when.limit
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
                        layout(mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
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
                           if (ed22.ci){
                              #------------------------------------------------------------#
                              #     Find the plot limit for the y scale.                   #
                              #------------------------------------------------------------#
                              yuse   = c(ed2.q025[d,],ed2.q975[d,],ed2.expected[d,]
                                        ,sta.q025[d,],sta.q975[d,],sta.expected[d,])
                              ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #    Make the polygons.                                      #
                              #------------------------------------------------------------#
                              size.poly     = list()
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
                              yuse   = c(ed2.expected[d,],sta.expected[d,]
                                        ,sta.q025[d,],sta.q975[d,])
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


                           #---------------------------------------------------------------#
                           #     Over-write with global rate.                              #
                           #---------------------------------------------------------------#
                           if (global.ylim){
                              ylt.mod = ylt[[ed2.rate]][[indiv[i]]]$global
                              ylt.min = mult * ylt.mod$min[d]
                              ylt.max = mult * ylt.mod$max[d]
                              ylimit  = pretty.xylim(u=c(ylt.min,ylt.max)
                                                    ,fracexp=0.0,is.log=ylog)
                           }#end if
                           #---------------------------------------------------------------#


                           #----- Set up the title and axes labels. -----------------------#
                           lesub = paste("DBH class:",dbh.names[d],sep="")
                           #---------------------------------------------------------------#


                           #----- Plot the box plot. --------------------------------------#
                           par(mar=lo.box$mar0)
                           #----- Plotting window and grid. -------------------------------#
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                           axis(side=2,las=1)
                           box()
                           title(main=lesub,xlab="",ylab="")
                           if (plotgrid){
                              abline(h=axTicks(2),v=when4,col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if
                           #----- Plot the taxon rate with confidence interval. -----------#
                           epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                                   ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
                           lines (x=whenmid4,y=sta.expected[d,],lwd=2.0,col=col.sta[1])
                           lines (x=whenmid4,y=ed2.expected[d,],lwd=2.0,col=col.ed2[1])
                           points(x=whenmid4,y=sta.expected[d,],pch=16 ,col=col.sta[1])
                           points(x=whenmid4,y=ed2.expected[d,],pch=16 ,col=col.ed2[1])
                           #---------------------------------------------------------------#
                        }#end for (d in 1:n.dbh)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(desc.rate," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=untab$cm)
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
                              , cex.main  = 0.9*cex.ptsz
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






                  #========================================================================#
                  #========================================================================#
                  #     PFT-dependent rates.                                               #
                  #------------------------------------------------------------------------#
                  if (pfttoo){
                     cat("  - PFT-dependent, plot level: ",desc.indiv[i],"...","\n")
                     ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
                     sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outrate   = paste(pftplot,ed2.rate,sep="/")
                     if (! file.exists(outrate)) dir.create(outrate)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outrate,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#


                     #----- Load the modelled rates. --------------------------------------#
                     mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                     sta.mod      = sta[[sta.rate]][[indiv[i]]]
                     sta.global   = sta[[sta.rate]][[indiv[i]]]$global
                     sta.expected = mult * rbind(sta.mod$expected[,-1],sta.global[1,-1])
                     sta.q025     = mult * rbind(sta.mod$q025    [,-1],sta.global[2,-1])
                     sta.q975     = mult * rbind(sta.mod$q975    [,-1],sta.global[3,-1])
                     ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]
                     ed2.global   = ed2[[ed2.rate]][[indiv[i]]]$global
                     ed2.expected = mult * rbind(ed2.mod$expected[,-1],ed2.global[1,-1])
                     ed2.q025     = mult * rbind(ed2.mod$q025    [,-1],ed2.global[2,-1])
                     ed2.q975     = mult * rbind(ed2.mod$q975    [,-1],ed2.global[3,-1])
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Find the DBH for x scale.                                        #
                     #---------------------------------------------------------------------#
                     x.years  = whenmid4
                     xlimit   = when.limit
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #    Find the limits for Y scale.                                     #
                     #---------------------------------------------------------------------#
                     if (ed22.ci){
                        #----- Use confidence intervals for limits. -----------------------#
                        yuse   = c(sta.q025,sta.q975,sta.expected
                                  ,ed2.q025,ed2.q975,ed2.expected)
                        ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)
                        #------------------------------------------------------------------#
                     }else{
                        #----- Use ED2 expected value for limits. -------------------------#
                        yuse   = c(sta.q025,sta.q975,sta.expected,ed2.expected)
                        ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#



                     #----- Find the appropriate number of boxes for PFTs. ----------------#
                     lo.pft = pretty.box(npfts+1,horizontal=TRUE)
                     #---------------------------------------------------------------------#




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



                        #------------------------------------------------------------------#
                        #     Split the window.                                            #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par.orig = par(no.readonly = TRUE)
                        par(oma = c(0.2,3,4,0))
                        layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                              ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#



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
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Loop over the PFTs.                                          #
                        #------------------------------------------------------------------#
                        for (p in sequence(npfts+1)){
                           if (p == npfts + 1){
                              pft.label="All PFTs"
                           }else{
                              pft.label=pft$name[mypfts[p]]
                           }#end if


                           #---------------------------------------------------------------#
                           #     Over-write with global rate.                              #
                           #---------------------------------------------------------------#
                           if (global.ylim){
                              ylt.mod = ylt[[ed2.rate]][[indiv[i]]]$taxon
                              ylt.min = mult * ylt.mod$min
                              ylt.max = mult * ylt.mod$max
                              ylimit  = pretty.xylim(u=c(ylt.min,ylt.max)
                                                    ,fracexp=0.0,is.log=ylog)
                              #------------------------------------------------------------#
                              #     Define margins.                                        #
                              #------------------------------------------------------------#
                              left    = TRUE
                              top     = TRUE
                              mar.now = lo.pft$mar0
                              #------------------------------------------------------------#
                           }else{


                              #------------------------------------------------------------#
                              #     Define margins.                                        #
                              #------------------------------------------------------------#
                              left    = lo.pft$left [p ]
                              top     = lo.pft$right[p ]
                              mar.now = lo.pft$mar  [p,]
                              #------------------------------------------------------------#
                           }#end if
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #    Make the polygons.                                         #
                           #---------------------------------------------------------------#
                           if (ed22.ci){
                              plot.poly     = list()
                              plot.poly$x   = c(x.years ,rev(x.years) ,NA
                                               ,x.years ,rev(x.years) )
                              plot.poly$y   = c(sta.q025[p,],rev(sta.q975[p,]),NA
                                               ,ed2.q025[p,],rev(ed2.q975[p,]))
                              plot.poly$col = c(col.sta[2],col.ed2[2])
                           }else{
                              plot.poly     = list()
                              plot.poly$x   = c(x.years ,rev(x.years) )
                              plot.poly$y   = c(sta.q025[p,],rev(sta.q975[p,]))
                              plot.poly$col = col.sta[2]
                           }#end if
                           #---------------------------------------------------------------#




                           #---------------------------------------------------------------#
                           #     Go on and plot stuff.                                     #
                           #---------------------------------------------------------------#
                           par(mar=mar.now)
                           #----- Plotting window and grid. -------------------------------#
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                           box()
                           axis.rt(side=1,at=when4,labels=when.label,las=5,off=.05)
                           axis(side=2,las=1)
                           title(main=pft.label)
                           if (plotgrid){
                              abline(h=axTicks(2),v=when4,col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if

                           #----- Plot the taxon rate with confidence interval. -----------#
                           epolygon(x=plot.poly$x,y=plot.poly$y,col=plot.poly$col
                                   ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
                           lines (x=x.years,y=sta.expected[p,],col=col.sta[1],lwd=2.0)
                           lines (x=x.years,y=ed2.expected[p,],col=col.ed2[1],lwd=2.0)
                           points(x=x.years,y=sta.expected[p,],col=col.sta[1],pch=16 )
                           points(x=x.years,y=ed2.expected[p,],col=col.ed2[1],pch=16 )
                           #---------------------------------------------------------------#
                        }#end for (p in sequence(npfts))
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(desc.rate," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
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
                              , cex.main  = 0.9*cex.ptsz
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
                     }#end for (o in 1:nout)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo)
                  #========================================================================#
                  #========================================================================#







                  #========================================================================#
                  #========================================================================#
                  #     PFT-dependent, DBH-LEVEL rates.                                    #
                  #------------------------------------------------------------------------#
                  if (pfttoo && sizetoo){
                     cat("  - PFT-dependent, DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outrate   = paste(pftsize,ed2.rate,sep="/")
                     if (! file.exists(outrate)) dir.create(outrate)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outrate,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#


                     #----- Load the modelled rates. --------------------------------------#
                     mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                     sta.mod      = sta[[sta.rate]][[indiv[i]]]
                     sta.expected = mult * sta.mod$expected[,,-1]
                     sta.q025     = mult * sta.mod$q025    [,,-1]
                     sta.q975     = mult * sta.mod$q975    [,,-1]
                     ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]
                     ed2.expected = mult * ed2.mod$expected[,,-1]
                     ed2.q025     = mult * ed2.mod$q025    [,,-1]
                     ed2.q975     = mult * ed2.mod$q975    [,,-1]
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Find the DBH for x scale.                                        #
                     #---------------------------------------------------------------------#
                     x.years  = whenmid4
                     xlimit   = when.limit
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Loop over PFTs.                                                 #
                     #---------------------------------------------------------------------#
                     for (p in sequence(npfts)){
                        pft.label = pft$name[mypfts[p]]
                        pft.key   = paste("pft",sprintf("%2.2i",mypfts[p]),sep="")


                        #------------------------------------------------------------------#
                        #     Loop over all formats, and make the plots.                   #
                        #------------------------------------------------------------------#
                        for (o in 1:nout){
                           #----- Open the file or the plot window. -----------------------#
                           fichier = file.path(outindiv
                                              ,paste(iata,"-tseries-",ed2.rate,"-",indiv[i]
                                                    ,"-",pft.key,"-",iphen.key[ph]
                                                    ,".",outform[o],sep=""))
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
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Split the window into several smaller windows.            #
                           #---------------------------------------------------------------#
                           par(par.user)
                           par.orig = par(no.readonly = TRUE)
                           par(oma = c(0.2,3,4,0))
                           layout(mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
                                 ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                                 )#end layout
                           #---------------------------------------------------------------#



                           #----- Plot legend. --------------------------------------------#
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
                           #---------------------------------------------------------------#




                           #---------------------------------------------------------------#
                           #    Loop over all DBH classes.                                 #
                           #---------------------------------------------------------------#
                           for (d in 1:n.dbh){

                              if (ed22.ci){
                                 #---------------------------------------------------------#
                                 #     Find the plot limit for the y scale.                #
                                 #---------------------------------------------------------#
                                 yuse   = c(ed2.q025[p,d,],ed2.q975[p,d,]
                                           ,sta.expected[p,d,]
                                           ,sta.q025[p,d,],sta.q975[p,d,]
                                           ,ed2.expected[p,d,])
                                 ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                                 #---------------------------------------------------------#



                                 #---------------------------------------------------------#
                                 #    Make the polygons.                                   #
                                 #---------------------------------------------------------#
                                 size.poly  = list()
                                 size.poly$x   = c(x.years     ,rev(x.years)     ,NA
                                                  ,x.years     ,rev(x.years)         )
                                 size.poly$y   = c(sta.q025[p,d,],rev(sta.q975[p,d,]),NA
                                                  ,ed2.q025[p,d,],rev(ed2.q975[p,d,])    )
                                 size.poly$col = c(col.sta[2]  ,col.ed2[2]           )
                                 #---------------------------------------------------------#
                              }else{
                                 #---------------------------------------------------------#
                                 #     Find the plot limit for the y scale.                #
                                 #---------------------------------------------------------#
                                 yuse   = c(ed2.expected[p,d,],sta.expected[d,]
                                           ,sta.q025[p,d,],sta.q975[p,d,])
                                 ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                                 #---------------------------------------------------------#



                                 #---------------------------------------------------------#
                                 #    Make the polygons.                                   #
                                 #---------------------------------------------------------#
                                 size.poly     = list()
                                 size.poly$x   = c(x.years       ,rev(x.years)           )
                                 size.poly$y   = c(sta.q025[p,d,],rev(sta.q975[p,d,])    )
                                 size.poly$col = c(col.sta[2]  )
                                 #---------------------------------------------------------#
                              }#end if
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #     Over-write with global rate.                           #
                              #------------------------------------------------------------#
                              if (global.ylim){
                                 ylt.mod = ylt[[ed2.rate]][[indiv[i]]]$taxon
                                 ylt.min = mult * ylt.mod$min[p,d]
                                 ylt.max = mult * ylt.mod$max[p,d]
                                 ylimit  = pretty.xylim(u=c(ylt.min,ylt.max)
                                                       ,fracexp=0.0,is.log=ylog)
                              }#end if
                              #------------------------------------------------------------#


                              #----- Set up the title and axes labels. --------------------#
                              lesub = paste("DBH class:",dbh.names[d],sep="")
                              #------------------------------------------------------------#


                              #----- Plot the box plot. -----------------------------------#
                              par(mar=lo.box$mar0)
                              #----- Plotting window and grid. ----------------------------#
                              plot.new()
                              plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                              axis(side=1)
                              axis(side=2,las=1)
                              box()
                              title(main=lesub,xlab="",ylab="")
                              if (plotgrid){
                                 grid(col=grid.colour,lwd=0.75,lty="dotted")
                              }#end if
                              #----- Plot the taxon rate with confidence interval. --------#
                              epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                                      ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
                              lines (x=x.years,y=sta.expected[p,d,],lwd=2.0,col=col.sta[1])
                              lines (x=x.years,y=ed2.expected[p,d,],lwd=2.0,col=col.ed2[1])
                              points(x=x.years,y=sta.expected[p,d,],pch=16 ,col=col.sta[1])
                              points(x=x.years,y=ed2.expected[p,d,],pch=16 ,col=col.ed2[1])
                              #------------------------------------------------------------#
                           }#end for (d in 1:n.dbh)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Make the title and axis labels.                           #
                           #---------------------------------------------------------------#
                           letitre = paste(desc.rate," - ",longname," - ",pft.label
                                          ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                           ley     = desc.unit(desc=desc.rate,unit=unit.rate[i])
                           lex     = "" # desc.unit(desc="Census year",unit=untab$cm)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Split the plotting window.                                #
                           #---------------------------------------------------------------#
                           gtitle( main      = letitre
                                 , xlab      = lex
                                 , ylab      = ley
                                 , off.xlab  = 1/6
                                 , line.xlab = 4.1
                                 , line.ylab = 2.6
                                 , cex.main  = 0.9*cex.ptsz
                                 )#end gtitle
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for (o in 1:nout)
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts))
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo && sizetoo)
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #==============================================================================#
            #==============================================================================#





            #==============================================================================#
            #==============================================================================#
            #  10.  Plot the theme time series of the rates (by DBH class if applicable).  #
            #==============================================================================#
            #==============================================================================#
            for (n in 1:npratethbpw){
               this.plot  = pratethbpw[[n]]
               sizetoo    = this.plot$sizetoo
               pfttoo     = this.plot$pfttoo
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
               col.ed2    = matrix(col.ed2,nrow=nrate)

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
                  outthbpw   = paste(outplot,theme.now,sep="/")
                  if (! file.exists(outthbpw)) dir.create(outthbpw)
                  #------------------------------------------------------------------------#


                  #----- Create path for this individual. ---------------------------------#
                  outindiv   = paste(outthbpw,indiv[i],sep="/")
                  if (! file.exists(outindiv)) dir.create(outindiv)
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #      Define the limits between boxes.                                  #
                  #------------------------------------------------------------------------#
                  off        = 0.1
                  width.bar  = (when4[-1] - when4[-n.census])
                  width.year = (nrate + 1) * width.bar + 1
                  x.bounds   = c(0, cumsum(width.year))
                  xlimit     = range(x.bounds)
                  #------------------------------------------------------------------------#


                  #----- Load the modelled rates. -----------------------------------------#
                  mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                  x.left       = list()
                  x.right      = list()
                  x.mid        = list()
                  mod.expected = list()
                  mod.q025     = list()
                  mod.q975     = list()
                  mod.colour   = list()
                  for (r in sequence(nrate+1)){
                     if (r == nrate + 1){
                        mod.now           = sta[[sta.rate]][[indiv[i]]]$global
                        mod.colour[[r]]   = col.sta
                     }else{
                        mod.now           = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                        mod.colour[[r]]   = col.ed2[r]
                     }#end if

                     x.left      [[r]] = ( x.bounds[-n.census]
                                         + 0.5 + (r-1) * width.bar + off)
                     x.right     [[r]] = x.left[[r]] + width.bar - 2. * off
                     x.mid       [[r]] = 0.5 * ( x.left[[r]] + x.right[[r]] )
                     mod.expected[[r]] = mult * mod.now[1,-1]
                     if (ed22.ci || r == nrate + 1){
                        mod.q025 [[r]] = mult * mod.now[2,-1]
                        mod.q975 [[r]] = mult * mod.now[3,-1]
                     }else{
                        mod.q025 [[r]] = NA   * mod.now[2,-1]
                        mod.q975 [[r]] = NA   * mod.now[3,-1]
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Find the plot limit for the y scale.                               #
                  #------------------------------------------------------------------------#
                  yuse   = c(0,unlist(mod.q025),unlist(mod.q975),unlist(mod.expected))
                  ylimit = pretty.xylim(u=yuse,fracexp=0.04,is.log=ylog)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Loop over all formats, and make the plots.                         #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Open the file or the plot window. -----------------------------#
                     fichier = file.path(outindiv
                                        ,paste(iata,"-thbpw-",theme.now,"-",indiv[i]
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
                     letitre = paste(theme.desc," - ",longname
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                     ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                     lex     = "" # desc.unit(desc="Census year",unit=NULL)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Go on and plot stuff.                                           #
                     #---------------------------------------------------------------------#
                     par(par.user)
                     layout(mat=rbind(2,1),heights=c(5,1))
                     #---------------------------------------------------------------------#


                     #----- Plot legend. --------------------------------------------------#
                     par(mar=c(0.1,0.1,0.1,0.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                     legend ( x       = "bottom"
                            , inset   = 0.0
                            , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                            , fill    = c(col.sta,rep("transparent",times=nrate-1)
                                         ,col.ed2)
                            , border  = c(foreground,rep("transparent",times=nrate-1)
                                         ,rep(foreground,times=nrate))
                            , density = -1
                            , ncol    = 2
                            , title   = "Whiskers - 95% C.I."
                            , cex     = 0.9 * cex.ptsz
                            , xpd     = TRUE
                            )#end legend
                     #---------------------------------------------------------------------#


                     #----- Plotting window and grid. -------------------------------------#
                     par(mar=c(5.1,4.4,4.1,2.1))
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,log=plog,yaxs="i")
                     axis.rt(side=1,at=x.bounds,labels=when.label,las=5,off=.05)
                     axis(side=2,las=1)
                     if (plotgrid){
                        abline(v=x.bounds,h=axTicks(2),col=grid.colour
                              ,lwd=0.75,lty="dotted")
                     }#end if
                     box()
                     title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.ptsz)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Plot rectangles and error bars.                                 #
                     #---------------------------------------------------------------------#
                     for (r in sequence(nrate+1)){
                        rect( xleft   = x.left[[r]]
                            , ybottom = 0
                            , xright  = x.right[[r]]
                            , ytop    = mod.expected[[r]]
                            , col     = mod.colour[[r]]
                            , border  = foreground
                            , lwd     = 1.5
                            )#end rect
                        if (ed22.ci || r == nrate+1){
                           whiskers( x       = x.mid[[r]]
                                   , y       = mod.expected[[r]]
                                   , ybottom = mod.q025[[r]]
                                   , ytop    = mod.q975[[r]]
                                   , col     = foreground
                                   , lwd     = 1.5
                                   )#end whiskers
                        }#end if
                        #------------------------------------------------------------------#
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
                  }#end for (o in 1:nout)
                  #------------------------------------------------------------------------#








                  #========================================================================#
                  #========================================================================#
                  #     DBH-LEVEL rates.                                                   #
                  #------------------------------------------------------------------------#
                  if (sizetoo){
                     cat("  - DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.part = paste(this.plot$ed2.part,"size",sep=".")
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
                     #      Define the limits between boxes.                               #
                     #---------------------------------------------------------------------#
                     off        = 0.1
                     width.bar  = (when4[-1] - when4[-n.census])
                     width.year = (nrate + 1) * width.bar + 1
                     x.bounds   = c(0, cumsum(width.year))
                     xlimit     = range(x.bounds)
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
                                           ,paste(iata,"-thbpw-",theme.now,"-",indiv[i],"-"
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
                        layout(mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
                              ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#


                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0.0
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta,rep("transparent",times=nrate-1)
                                            ,col.ed2)
                               , border  = c(foreground,rep("transparent",times=nrate-1)
                                            ,rep(foreground,times=nrate))
                               , density = -1
                               , ncol    = 2
                               , title   = "Whiskers - 95% C.I."
                               , cex     = cex.ptsz
                               , xpd     = TRUE
                               )#end legend
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #    Loop over all DBH classes.                                    #
                        #------------------------------------------------------------------#
                        for (d in 1:n.dbh){



                           #----- Load the modelled rates. --------------------------------#
                           mult = 100 - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                           x.left       = list()
                           x.right      = list()
                           x.mid        = list()
                           mod.expected = list()
                           mod.q025     = list()
                           mod.q975     = list()
                           mod.colour   = list()
                           for (r in sequence(nrate+1)){
                              if (r == nrate + 1){
                                 mod.now           = sta[[sta.rate]][[indiv[i]]]$global
                                 mod.colour[[r]]   = col.sta
                              }else{
                                 mod.now           = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                                 mod.colour[[r]]   = col.ed2[r]
                              }#end if

                              x.left      [[r]] = ( x.bounds[-n.census]
                                                  + 0.5 + (r-1) * width.bar + off)
                              x.right     [[r]] = x.left[[r]] + width.bar - 2. * off
                              x.mid       [[r]] = 0.5 * ( x.left[[r]] + x.right[[r]] )
                              mod.expected[[r]] = mult * mod.now[1,d,-1]
                              if (ed22.ci || r == nrate + 1){
                                 mod.q025 [[r]] = mult * mod.now[2,d,-1]
                                 mod.q975 [[r]] = mult * mod.now[3,d,-1]
                              }else{
                                 mod.q025 [[r]] = NA   * mod.now[2,d,-1]
                                 mod.q975 [[r]] = NA   * mod.now[3,d,-1]
                              }#end if
                              #------------------------------------------------------------#
                           }#end for
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Find the plot limit for the y scale.                      #
                           #---------------------------------------------------------------#
                           yuse   = c(0,unlist(mod.q025),unlist(mod.q975)
                                     ,unlist(mod.expected))
                           ylimit = pretty.xylim(u=yuse,fracexp=0.04,is.log=ylog)
                           #---------------------------------------------------------------#


                           #----- Set up the title and axes labels. -----------------------#
                           lesub = paste("DBH class:",dbh.names[d],sep="")
                           #---------------------------------------------------------------#


                           #----- Plotting window and grid. -------------------------------#
                           par(mar=lo.box$mar0)
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog,yaxs="i")
                           axis.rt(side=1,at=x.bounds,labels=when.label
                                  ,las=5,off=.05)
                           axis(side=2,las=1)
                           if (plotgrid){
                              abline(v=x.bounds,h=axTicks(2),col=grid.colour,lty="dotted")
                           }#end if
                           box()
                           title(main=lesub,xlab="",ylab="")
                           #---------------------------------------------------------------#


                           #---------------------------------------------------------------#
                           #     Plot rectangles and error bars.                           #
                           #---------------------------------------------------------------#
                           for (r in sequence(nrate+1)){
                              rect( xleft   = x.left[[r]]
                                  , ybottom = 0
                                  , xright  = x.right[[r]]
                                  , ytop    = mod.expected[[r]]
                                  , col     = mod.colour[[r]]
                                  , border  = foreground
                                  , lwd     = 1.5
                                  )#end rect
                              if (ed22.ci || r == nrate+1){
                                 whiskers( x       = x.mid[[r]]
                                         , y       = mod.expected[[r]]
                                         , ybottom = mod.q025[[r]]
                                         , ytop    = mod.q975[[r]]
                                         , col     = foreground
                                         , lwd     = 1.5
                                         )#end whiskers
                              }#end if
                              #------------------------------------------------------------#
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (d in 1:n.dbh)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
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
                              , cex.main  = 0.9*cex.ptsz
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





                  #========================================================================#
                  #========================================================================#
                  #    PFT-specific plots.                                                 #
                  #------------------------------------------------------------------------#
                  if (pfttoo){

                     #=====================================================================#
                     #=====================================================================#
                     #     PLOT-LEVEL rates.                                               #
                     #---------------------------------------------------------------------#
                     cat("  - PFT-dependent plot level: ",desc.indiv[i],"...","\n")
                     ed2.rate = paste(this.plot$ed2.rate,"plot",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"plot",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(pftplot,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #      Define the limits between boxes.                               #
                     #---------------------------------------------------------------------#
                     off        = 0.1
                     width.bar  = (when4[-1] - when4[-n.census])
                     width.year = (nrate + 1) * width.bar + 1
                     x.bounds   = c(0, cumsum(width.year))
                     xlimit     = range(x.bounds)
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.pft = pretty.box(n=npfts+1,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #----- Load the modelled rates. --------------------------------------#
                     mult         = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                     x.left       = list()
                     x.right      = list()
                     x.mid        = list()
                     mod.expected = list()
                     mod.q025     = list()
                     mod.q975     = list()
                     mod.colour   = list()
                     for (r in sequence(nrate+1)){
                        if (r == nrate + 1){
                           mod.now           = sta[[sta.rate]][[indiv[i]]]
                           mod.global        = sta[[sta.rate]][[indiv[i]]]$global
                           mod.colour[[r]]   = col.sta
                        }else{
                           mod.now           = ed2[[ed2.rate[r]]][[indiv[i]]]
                           mod.global        = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                           mod.colour[[r]]   = col.ed2[r]
                        }#end if

                        x.left      [[r]] = ( x.bounds[-n.census]
                                            + 0.5 + (r-1) * width.bar + off)
                        x.right     [[r]] = x.left[[r]] + width.bar - 2. * off
                        x.mid       [[r]] = 0.5 * ( x.left[[r]] + x.right[[r]] )
                        mod.expected[[r]] = mult * rbind(mod.now$expected[ ,-1]
                                                        ,mod.global      [1,-1])
                        if (ed22.ci || r == nrate + 1){
                           mod.q025    [[r]] = mult * rbind(mod.now$q025    [ ,-1]
                                                           ,mod.global      [2,-1])
                           mod.q975    [[r]] = mult * rbind(mod.now$q975    [ ,-1]
                                                           ,mod.global      [3,-1])
                        }else{
                           mod.q025    [[r]] = NA   * rbind(mod.now$q025    [ ,-1]
                                                           ,mod.global      [2,-1])
                           mod.q975    [[r]] = NA   * rbind(mod.now$q975    [ ,-1]
                                                           ,mod.global      [3,-1])
                        }#end if
                        #---------------------------------------------------------------------#
                     }#end for
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Find Y limits.                                                  #
                     #---------------------------------------------------------------------#
                     yuse   = c(0,unlist(mod.q025),unlist(mod.q975),unlist(mod.expected))
                     ylimit = pretty.xylim(u=yuse,fracexp=0.04,is.log=ylog)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Loop over all formats, and make the plots.                      #
                     #---------------------------------------------------------------------#
                     for (o in 1:nout){
                        #----- Open the file or the plot window. --------------------------#
                        fichier = file.path(outindiv
                                           ,paste(iata,"-thbpw-",theme.now,"-",indiv[i]
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
                           pdf(file=fichier,onefile=FALSE,width=size$width
                              ,height=size$height,pointsize=ptsz,paper=size$paper)
                        }#end if
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Split the window.                                            #
                        #------------------------------------------------------------------#
                        par(par.user)
                        par.orig = par(no.readonly = TRUE)
                        par(oma = c(0.2,3,4,0))
                        layout(mat    = rbind(lo.pft$mat.off,rep(1,times=lo.pft$ncol))
                              ,height = c(rep(5/lo.pft$nrow,times=lo.pft$nrow),1)
                              )#end layout
                        #------------------------------------------------------------------#


                        #----- Plot legend. -----------------------------------------------#
                        par(mar=c(0.1,0.1,0.1,0.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                        legend ( x       = "bottom"
                               , inset   = 0.0
                               , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                               , fill    = c(col.sta,rep("transparent",times=nrate-1)
                                            ,col.ed2)
                               , border  = c(foreground,rep("transparent",times=nrate-1)
                                            ,rep(foreground,times=nrate))
                               , density = -1
                               , ncol    = 2
                               , title   = "Whiskers - 95% C.I."
                               , cex     = cex.ptsz
                               , xpd     = TRUE
                               )#end legend
                        #------------------------------------------------------------------#




                        #------------------------------------------------------------------#
                        #     Loop over all PFTS to be plotted.                            #
                        #------------------------------------------------------------------#
                        for (p in sequence(npfts+1)){
                           if (p == npfts + 1){
                              pft.label = "All PFTs"
                           }else{
                              pft.label = pft$name[mypfts[p]]
                           }#end if


                           #---------------------------------------------------------------#
                           #     Define margins.                                           #
                           #---------------------------------------------------------------#
                           left    = lo.pft$left [p ]
                           right   = lo.pft$right[p ]
                           mar.now = lo.pft$mar  [p,]
                           #---------------------------------------------------------------#




                           #----- Plotting window and grid. -------------------------------#
                           par(mar=c(4.1,3.1,2.1,2.1))
                           plot.new()
                           plot.window(xlim=xlimit,ylim=ylimit,log=plog,yaxs="i")
                           box()
                           axis.rt(side=1,at=x.bounds,labels=when.label
                                  ,las=5,off=.05)
                           axis(side=2,las=1)
                           if (plotgrid){
                              abline(v=x.bounds,h=axTicks(2),col=grid.colour
                                    ,lwd=0.75,lty="dotted")
                           }#end if
                           title(main=pft.label)
                           #---------------------------------------------------------------#




                           #---------------------------------------------------------------#
                           #     Plot rectangles and error bars.                           #
                           #---------------------------------------------------------------#
                           for (r in sequence(nrate+1)){
                              rect( xleft   = x.left[[r]]
                                  , ybottom = 0.
                                  , xright  = x.right[[r]]
                                  , ytop    = mod.expected[[r]][p,]
                                  , col     = mod.colour[[r]]
                                  , border  = foreground
                                  , lwd     = 1.5
                                  )
                              if (ed22.ci || r == nrate+1){
                                 whiskers( x       = x.mid[[r]]
                                         , y       = mod.expected[[r]][p,]
                                         , ybottom = mod.q025[[r]][p,]
                                         , ytop    = mod.q975[[r]][p,]
                                         , col     = foreground
                                         , lwd     = 1.5
                                         )
                              }#end if
                              #------------------------------------------------------------#
                           }#end for
                           #---------------------------------------------------------------#
                        }#end for (p in sequence(npfts)
                        #------------------------------------------------------------------#






                        #------------------------------------------------------------------#
                        #     Make the title and axis labels.                              #
                        #------------------------------------------------------------------#
                        letitre = paste(theme.desc," - ",longname
                                       ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                        ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                        lex     = "" # desc.unit(desc="Census year",unit=NULL)
                        #------------------------------------------------------------------#


                        #------------------------------------------------------------------#
                        #      Plot global annotation.                                     #
                        #------------------------------------------------------------------#
                        gtitle( main      = letitre
                              , xlab      = lex
                              , ylab      = ley
                              , off.xlab  = 1/6
                              , line.xlab = 4.1
                              , line.ylab = 2.6
                              , cex.main  = 0.9*cex.ptsz
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
                     }#end for (o in 1:nout)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo)
                  #========================================================================#
                  #========================================================================#





                  #========================================================================#
                  #========================================================================#
                  #     DBH-LEVEL rates by PFT.                                            #
                  #------------------------------------------------------------------------#
                  if (sizetoo && pfttoo){
                     cat("  - PFT-dependent, DBH classes: ",desc.indiv[i],"...","\n")
                     ed2.part = paste(this.plot$ed2.part,"size",sep=".")
                     ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
                     sta.rate = paste(this.plot$sta.rate,"size",sep=".")


                     #----- Create a directory for this type of plot. ---------------------#
                     outtheme   = paste(pftsize,theme.now,sep="/")
                     if (! file.exists(outtheme)) dir.create(outtheme)
                     #---------------------------------------------------------------------#


                     #----- Create path for this individual. ------------------------------#
                     outindiv   = paste(outtheme,indiv[i],sep="/")
                     if (! file.exists(outindiv)) dir.create(outindiv)
                     #---------------------------------------------------------------------#





                     #---------------------------------------------------------------------#
                     #      Define the limits between boxes.                               #
                     #---------------------------------------------------------------------#
                     off        = 0.1
                     width.bar  = (when4[-1] - when4[-n.census])
                     width.year = (nrate + 1) * width.bar + 1
                     x.bounds   = c(0, cumsum(width.year))
                     xlimit     = range(x.bounds)
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #      Define a nice configuration for the multiple panels.           #
                     #---------------------------------------------------------------------#
                     lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #    Loop over all PFTs.                                              #
                     #---------------------------------------------------------------------#
                     for (p in sequence(npfts)){
                        pft.label = pft$name[mypfts[p]]
                        pft.key   = paste("pft",sprintf("%2.2i",mypfts[p]),sep="")



                        #------------------------------------------------------------------#
                        #     Loop over all formats, and make the plots.                   #
                        #------------------------------------------------------------------#
                        for (o in 1:nout){
                           #----- Open the file or the plot window. -----------------------#
                           fichier = file.path(outindiv
                                              ,paste(iata,"-thbpw-",theme.now,"-",indiv[i],
                                                    "-",pft.key,"-",iphen.key[ph]
                                                    ,".",outform[o] ,sep=""))
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
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Split the window into several smaller windows.            #
                           #---------------------------------------------------------------#
                           par(par.user)
                           par.orig = par(no.readonly = TRUE)
                           par(oma = c(0.2,3,4,0))
                           layout(mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
                                 ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                                 )#end layout
                           #---------------------------------------------------------------#


                           #----- Plot legend. --------------------------------------------#
                           par(mar=c(0.1,0.1,0.1,0.1))
                           plot.new()
                           plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                           legend ( x       = "bottom"
                                  , inset   = 0.0
                                  , legend  = c("Observed",rep("",times=nrate-1),desc.rate)
                                  , fill    = c(col.sta,rep("transparent",times=nrate-1)
                                               ,col.ed2)
                                  , border  = c(foreground,rep("transparent",times=nrate-1)
                                               ,rep(foreground,times=nrate))
                                  , density = -1
                                  , ncol    = 2
                                  , title   = "Whiskers - 95% C.I."
                                  , cex     = cex.ptsz
                                  , xpd     = TRUE
                                  )#end legend
                           #---------------------------------------------------------------#




                           #---------------------------------------------------------------#
                           #    Loop over all DBH classes.                                 #
                           #---------------------------------------------------------------#
                           for (d in 1:n.dbh){

                              #----- Load the modelled rates. -----------------------------#
                              mult = 100. - 99 * as.numeric(indiv[i] %in% c("acc","anpp"))
                              x.left       = list()
                              x.right      = list()
                              x.mid        = list()
                              mod.expected = list()
                              mod.q025     = list()
                              mod.q975     = list()
                              mod.colour   = list()
                              for (r in sequence(nrate+1)){
                                 if (r == nrate + 1){
                                    mod.now         = sta[[sta.rate]][[indiv[i]]]
                                    mod.colour[[r]] = col.sta
                                 }else{
                                    mod.now         = ed2[[ed2.rate[r]]][[indiv[i]]]
                                    mod.colour[[r]] = col.ed2[r]
                                 }#end if

                                 x.left      [[r]] = ( x.bounds[-n.census]
                                                     + 0.5 + (r-1) * width.bar + off)
                                 x.right     [[r]] = x.left[[r]] + width.bar - 2. * off
                                 x.mid       [[r]] = 0.5 * ( x.left[[r]] + x.right[[r]] )
                                 mod.expected[[r]] = mult * mod.now$expected[p,d,-1]
                                 if (ed22.ci || r == nrate + 1){
                                    mod.q025    [[r]] = mult * mod.now$q025    [p,d,-1]
                                    mod.q975    [[r]] = mult * mod.now$q975    [p,d,-1]
                                 }else{
                                    mod.q025    [[r]] = NA   * mod.now$q025    [p,d,-1]
                                    mod.q975    [[r]] = NA   * mod.now$q975    [p,d,-1]
                                 }#end if
                                 #---------------------------------------------------------#
                              }#end for
                              #------------------------------------------------------------#


                              #------------------------------------------------------------#
                              #     Find Y limits.                                         #
                              #------------------------------------------------------------#
                              yuse   = c(0,unlist(mod.q025),unlist(mod.q975)
                                          ,unlist(mod.expected))
                              ylimit = pretty.xylim(u=yuse,fracexp=0.04,is.log=ylog)
                              #------------------------------------------------------------#


                              #----- Set up the title and axes labels. --------------------#
                              lesub = paste("DBH class:",dbh.names[d],sep="")
                              #------------------------------------------------------------#


                              #----- Plotting window and grid. ----------------------------#
                              par(mar=lo.box$mar0)
                              plot.new()
                              plot.window(xlim=xlimit,ylim=ylimit,log=plog,yaxs="i")
                              axis.rt(side=1,at=x.bounds,labels=when.label
                                     ,las=5,off=.05)
                              axis(side=2,las=1)
                              if (plotgrid){
                                 abline(v=x.bounds,h=axTicks(2),col=grid.colour
                                       ,lwd=0.75,lty="dotted")
                              }#end if
                              box()
                              title(main=lesub,xlab="",ylab="")
                              #------------------------------------------------------------#



                              #------------------------------------------------------------#
                              #     Plot rectangles and error bars.                        #
                              #------------------------------------------------------------#
                              for (r in sequence(nrate+1)){
                                 rect( xleft   = x.left      [[r]]
                                     , ybottom = 0
                                     , xright  = x.right     [[r]]
                                     , ytop    = mod.expected[[r]]
                                     , col     = mod.colour  [[r]]
                                     , border  = foreground
                                     , lwd     = 1.5
                                     )#end rect
                                 if (ed22.ci || r == nrate+1){
                                    whiskers( x       = x.mid       [[r]]
                                            , y       = mod.expected[[r]]
                                            , ybottom = mod.q025    [[r]]
                                            , ytop    = mod.q975    [[r]]
                                            , col     = foreground
                                            , lwd     = 1.5
                                            )
                                 }#end if (ed22.ci || r == nrate+1)
                                 #---------------------------------------------------------#
                              }#end for (r in sequence(nrate+1))
                              #------------------------------------------------------------#
                           }#end for (d in 1:n.dbh)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Make the title and axis labels.                           #
                           #---------------------------------------------------------------#
                           letitre = paste(theme.desc," - ",longname," - ",pft.label
                                    ,"\n",iphen.desc[ph]," - ",tfall.desc[tf],sep="")
                           ley     = desc.unit(desc=theme.desc,unit=unit.rate[i])
                           lex     = "" # desc.unit(desc="Census year",unit=NULL)
                           #---------------------------------------------------------------#



                           #---------------------------------------------------------------#
                           #     Plot title.                                               #
                           #---------------------------------------------------------------#
                           gtitle( main      = letitre
                                 , xlab      = lex
                                 , ylab      = ley
                                 , off.xlab  = 1/6
                                 , line.xlab = 4.1
                                 , line.ylab = 2.6
                                 , cex.main  = 0.9*cex.ptsz
                                 )#end gtitle
                           #---------------------------------------------------------------#



                           #----- Close the device. ---------------------------------------#
                           if (outform[o] == "x11"){
                              locator(n=1)
                              dev.off()
                           }else{
                              dev.off()
                           }#end if
                           #---------------------------------------------------------------#
                        }#end for (o in 1:nout)
                        #------------------------------------------------------------------#
                     }#end for (p in sequence(npfts)
                     #---------------------------------------------------------------------#
                  }#end if (pfttoo && sizetoo)
                  #------------------------------------------------------------------------#

               }#end for (i in 1:nindiv)
               #---------------------------------------------------------------------------#
            }#end for (n in 1:npratethbpw)
            #==============================================================================#
            #==============================================================================#
         }#end if (census.name %in% ls())
         #=================================================================================#
         #=================================================================================#
      }#end for phenology
      #------------------------------------------------------------------------------------#

   }#end for tfall
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#
