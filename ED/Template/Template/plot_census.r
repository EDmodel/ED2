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
here           = "thispath"                          # Current directory.
there          = "thatpath"                          # Directory where analyses/history are 
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc" # Source  directory.
outroot        = "thisoutroot"                       # Directory for figures
#------------------------------------------------------------------------------------------#




#----- Time options. ----------------------------------------------------------------------#
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
biocyca        = mybiocyca         # First year of the met cycle
biocycz        = mybiocycz         # Last year of the met cycle
various.cycles = myvarcycle
sasmonth.short = c(2,5,8,11)
sasmonth.long  = 5
nyears.long    = 25
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("thispoly")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = thisoutform            # Formats for output file.  Supported formats are:
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
ncolsfc        = 200                    # Target number of colours for filled contour.
fcgrid         = TRUE                   # Include a grid on the filled contour plots?
ncolshov       = 200                    # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE                   # Include a grid on the Hovmoller plots?
legwhere       = "topleft"              # Where should I place the legend?
inset          = 0.01                   # Inset between legend and edge of plot region.
scalleg        = 0.40                   # Expand y limits by this relative amount to fit
                                        #    the legend
cex.main       = 0.8                    # Scale coefficient for the title
theta          = 315.                   # Azimuth for perspective projection
phi            = 30.                    # Vertical angle for perspective projection
ltheta         = -210.                  # Azimuth angle for light
shade          = 0.125                  # Shade intensity
expz           = 0.5                    # Expansion factor for Z axis
cexmin         = 0.5                    # Minimum "head" size of the lollipop
cexmax         = 3.0                    # Maximum "head" size of the lollipop
ylnudge         = 0.05                  # Nudging factor for ylimit
ptype          = "l"                    # Type of plot
ptyped         = "p"                    # Type of plot
ptypeb         = "o"                    # Type of plot
mtext.xoff     = -7.00                  # Offset for the x label
mtext.yoff     = -1.00                  # Offset for the y label
mtext.xadj     =  0.50                  # Offset for the x label
mtext.yadj     =  0.65                  # Offset for the y label
drought.mark   = mydroughtmark          # Put a background to highlight droughts?
drought.yeara  = mydroughtyeara         # First year that has drought
drought.yearz  = mydroughtyearz         # Last year that has drought
months.drought = mymonthsdrought        # Months with drought
ibackground    = mybackground           # Background colour
#------------------------------------------------------------------------------------------#


#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      = myidbhtype   # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
ed22.ci        = TRUE         # Plot confidence interval for ED?
n.boot         = 1000         # Number of realisations for bootstrap
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
recr.vars     = c("n")
recr.labels   = c("Individuals")
mort.vars     = c("n")
mort.labels   = c("Individuals")
growth.vars   = c("dbh","agb","ba")
growth.labels = c("DBH","Above Ground Biomass","Basal Area")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Comparisons.                                                                         #
#------------------------------------------------------------------------------------------#
#---- 1. Plot time series of median and confidence intervals. -----------------------------#
pratets      = list()
pratets[[1]] = list( ed2.rate   = "recr"
                   , sta.rate   = "recr"
                   , sizetoo    = FALSE
                   , pfttoo     = TRUE
                   , desc.rate  = "Recruitment rate"
                   , unit.rate  = "%/yr"
                   , col.ed2    = c(green.fg,green.bg)
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
                   , unit.rate  = "%/yr"
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
                   , unit.rate  = "%/yr"
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
                   , unit.rate  = "%/yr"
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
                   , unit.rate  = "%/yr"
                   , col.ed2    = c(yellow.fg,yellow.bg)
                   , col.sta    = c(grey.fg,grey.bg)
                   , indiv      = growth.vars
                   , desc.indiv = growth.labels
                   , legpos     = "topright"
                   , plog       = ""
                   )#end list



#---- 2. Plot median and confidence intervals for all size classes and censuses. ----------#
pratesize      = list()
pratesize[[1]] = list( ed2.rate   = "mort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Mortality rate"
                     , unit.rate  = "%/yr"
                     , col.ed2    = c(purple.fg,purple.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = mort.vars
                     , desc.indiv = mort.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[2]] = list( ed2.rate   = "ddmort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Density-dependent mort. rate"
                     , unit.rate  = "%/yr"
                     , col.ed2    = c(indigo.fg,indigo.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = mort.vars
                     , desc.indiv = mort.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[3]] = list( ed2.rate   = "dimort"
                     , sta.rate   = "mort"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Density-independent mort. rate"
                     , unit.rate  = "%/yr"
                     , col.ed2    = c(blue.fg,blue.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = mort.vars
                     , desc.indiv = mort.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
pratesize[[4]] = list( ed2.rate   = "growth"
                     , sta.rate   = "growth"
                     , sizetoo    = TRUE
                     , pfttoo     = TRUE
                     , desc.rate  = "Growth rate"
                     , unit.rate  = "%/yr"
                     , col.ed2    = c(yellow.fg,yellow.bg)
                     , col.sta    = c(grey.fg,grey.bg)
                     , indiv      = growth.vars
                     , desc.indiv = growth.labels
                     , legpos     = "topright"
                     , plog       = ""
                     )#end list
#---- 3. Plot median and confidence intervals for themes. ---------------------------------#
pratetheme      = list()
pratetheme[[1]] = list( ed2.rate   = c("ddmort","dimort")
                      , sta.rate   = "mort"
                      , sizetoo    = TRUE
                      , pfttoo     = TRUE
                      , desc.rate  = c("ED-2.2 Negative C Balance"
                                      ,"ED-2.2 Other")
                      , unit.rate  = "%/yr"
                      , col.ed2    = rbind( c(indigo.fg,indigo.bg)
                                          , c(green.fg ,green.bg)
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






#----- XYZ plots, to explore the parameter space. -----------------------------------------#
xyzvar = list()
xyzvar$xvar      = list( list( vname = "lai"
                             , desc  = "Leaf area index"
                             , unit  = "m2/m2"
                             , add   = 0.
                             , mult  = 1.
                             , leg   = "right"
                             , log   = FALSE
                             )#end list
                       , list( vname = "ba"
                             , desc  = "Basal area"
                             , unit  = "cm2/m2"
                             , add   = 0.
                             , mult  = 1.
                             , leg   = "right"
                             , log   = FALSE
                             )#end list
                       , list( vname = "agb"
                             , desc  = "Above-ground biomass"
                             , unit  = "kgC/m2"
                             , add   = 0.
                             , mult  = 1.
                             , leg   = "right"
                             , log   = FALSE
                             )#end list
                       )#end list
xyzvar$yvar      = list( list( vname      = "recr"
                             , desc       = "Recruitment rate"
                             , key        = "Recruitment"
                             , unit       = "%pop/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = FALSE
                             )#end list
                       , list( vname      = "mort"
                             , desc       = "Mortality rate"
                             , key        = "Mortality"
                             , unit       = "%pop/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = TRUE
                             )#end list
                       , list( vname      = "ddmort"
                             , desc       = "Density-dependent mortality rate"
                             , key        = "DD Mort."
                             , unit       = "%pop/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = TRUE
                             )#end list
                       , list( vname      = "dimort"
                             , desc       = "Density-independent mortality rate"
                             , key        = "DI Mort."
                             , unit       = "%pop/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = TRUE
                             )#end list
                       , list( vname      = "growdbh"
                             , desc       = "Growth rate (DBH)"
                             , key        = "Growth"
                             , unit       = "%DBH/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = TRUE
                             )#end list
                       , list( vname      = "growagb"
                             , desc       = "Growth rate (AGB)"
                             , key        = "Growth"
                             , unit       = "%AGB/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = TRUE
                             )#end list
                       , list( vname      = "growba"
                             , desc       = "Growth rate (BA)"
                             , key        = "Growth"
                             , unit       = "%BA/yr"
                             , add        = 0.
                             , mult       = 100.
                             , leg        = "top"
                             , log        = TRUE
                             , sizetoo    = TRUE
                             )#end list
                       )#end list
xyzvar$zvar      = list( list( vname      = "rshort"
                             , desc       = "Mean shortwave radiation"
                             , unit       = "W/m2"
                             , add        = 0.
                             , mult       = 1.
                             , col.scheme = "muitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "fs.open"
                             , desc       = "Minimum water stress scale"
                             , unit       = "--"
                             , add        = 0.
                             , mult       = 1.
                             , col.scheme = "imuitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "paw"
                             , desc       = "Minimum available water"
                             , unit       = "%"
                             , add        = 0.
                             , mult       = 100.
                             , col.scheme = "imuitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "smpot"
                             , desc       = "Maximum matric potential"
                             , unit       = "MPa"
                             , add        = 0.
                             , mult       = 1.
                             , col.scheme = "muitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "atm.vpd"
                             , desc       = "Above-canopy VPD"
                             , unit       = "hPa"
                             , add        = 0.
                             , mult       = 0.01
                             , col.scheme = "muitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "leaf.vpd"
                             , desc       = "Leaf-level VPD"
                             , unit       = "hPa"
                             , add        = 0.
                             , mult       = 0.01
                             , col.scheme = "muitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "leaf.gsw"
                             , desc       = "Stomatal conductance"
                             , unit       = "kg/m2/day"
                             , add        = 0.
                             , mult       = 86400.
                             , col.scheme = "imuitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "leaf.gbw"
                             , desc       = "Leaf Bnd. Lyr. Conductance"
                             , unit       = "kg/m2/day"
                             , add        = 0.
                             , mult       = 86400.
                             , col.scheme = "imuitas"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "cba.light"
                             , desc       = "Carbon balance (Max. light)"
                             , unit       = "kgC/m2/yr"
                             , add        = 0.
                             , mult       = 1.
                             , col.scheme = "clife"
                             , log        = FALSE
                             )#end list
                       , list( vname      = "cba.moist"
                             , desc       = "Carbon balance (Max. Moisture)"
                             , unit       = "kgC/m2/yr"
                             , add        = 0.
                             , mult       = 1.
                             , col.scheme = "clife"
                             , log        = FALSE
                             )#end list
                       )#end if
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


#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                      ,monthbeg=monthbeg)
   inpref  = thispoi$pathin
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"census",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz
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
   #     Reset the soil flag.                                                              #
   #---------------------------------------------------------------------------------------#
   read.soil = TRUE
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     We only run this part of the code if there are observations to compare with the   #
   # model.                                                                                #
   #---------------------------------------------------------------------------------------#
   if (census.name %in% ls()){

      #----- Check that the directories exist. --------------------------------------------#
      if (! file.exists(outmain)) dir.create(outmain)
      if (! file.exists(outpref)) dir.create(outpref)
      #------------------------------------------------------------------------------------#


      #----- Initialise the parameter space structure. ------------------------------------#
      pspace        = list()
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
      dbh.names   = dimnames(sta$mort.size$n$median)[[2]]
      year4       = numyears(sta$when)
      biocyca     = year4[2]
      biocycz     = year4[length(year4)]
      dyear       = c(NA,diff(year4))
      census.desc = paste(year4-c(NA,diff(year4)),year4,sep="-")


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
      #     Find out how many pseudo-census cycles we can repeat.  Notice that it is your  #
      # responsibility to make sure that met driver forcing and census cycles are          #
      # synchronised, this script will not check it for you.                               #
      #------------------------------------------------------------------------------------#
      act.census.yeara    = numyears(sta$when[2])
      act.census.yearz    = numyears(sta$when[n.census])
      act.census.year     = numyears(sta$when)
      #------------------------------------------------------------------------------------#

      ncyc                = biocycz - biocyca + 1
      if (various.cycles){
         nfullcyc         = 1 + floor( (act.census.yeara - 1 - yeara ) / ncyc )
      }else{
         nfullcyc         = 1
      }#end if
      i1stfull            = act.census.yeara - (nfullcyc - 1) * ncyc
      n.months            = length(census.idx)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Make the vector wit the environmental, plot-level, and size-level dimensions.   #
      #------------------------------------------------------------------------------------#
      dim.envr     = c(             n.months,nfullcyc)
      dim.plot     = c(npft+1,      n.months,nfullcyc)
      dim.size     = c(npft+1,n.dbh,n.months,nfullcyc)
      #------------------------------------------------------------------------------------#


      #----- Initialise all the structures for which we will compare. ---------------------#
      ts.n.plot    = array( 0., dim = dim.plot)
      ts.agb.plot  = array( NA, dim = dim.plot)
      ts.ba.plot   = array( NA, dim = dim.plot)
      ts.lai.plot  = array( NA, dim = dim.plot)
      ts.n.size    = array( 0., dim = dim.size)
      ts.agb.size  = array( NA, dim = dim.size)
      ts.ba.size   = array( NA, dim = dim.size)
      ts.lai.size  = array( NA, dim = dim.size)
      #----- Recruitment. -----------------------------------------------------------------#
      ts.recr.plot = list()
      for (v in 1:nrecr.vars){
         ts.recr.plot  [[recr.vars[v]]] = array( NA, dim = dim.plot)
      }#end for
      #----- Mortality. -------------------------------------------------------------------#
      ts.mort.plot     = list()
      ts.ddmort.plot   = list()
      ts.dimort.plot   = list()
      ts.mort.size     = list()
      ts.ddmort.size   = list()
      ts.dimort.size   = list()
      for (v in 1:nmort.vars){
         ts.mort.plot  [[mort.vars[v]]] = array( NA, dim = dim.plot)
         ts.ddmort.plot[[mort.vars[v]]] = array( NA, dim = dim.plot)
         ts.dimort.plot[[mort.vars[v]]] = array( NA, dim = dim.plot)
         ts.mort.size  [[mort.vars[v]]] = array( NA, dim = dim.size)
         ts.ddmort.size[[mort.vars[v]]] = array( NA, dim = dim.size)
         ts.dimort.size[[mort.vars[v]]] = array( NA, dim = dim.size)
      }#end for
      #----- Growth. ----------------------------------------------------------------------#
      ts.growth.plot   = list()
      ts.growth.size   = list()
      for (v in 1:ngrowth.vars){
         ts.growth.plot[[growth.vars[v]]] = array( NA, dim = dim.plot)
         ts.growth.size[[growth.vars[v]]] = array( NA, dim = dim.size)
      }#end for
      #----- Environmental variables. -----------------------------------------------------#
      ts.rshort        = array( NA, dim = dim.envr)
      ts.fs.open       = array( NA, dim = dim.envr)
      ts.paw           = array( NA, dim = dim.envr)
      ts.smpot         = array( NA, dim = dim.envr)
      ts.census.year   = array( NA, dim = dim.envr)
      ts.census.idx    = array( NA, dim = dim.envr)
      ts.atm.vpd       = array( NA, dim = dim.envr)
      ts.leaf.vpd      = array( NA, dim = dim.envr)
      ts.leaf.gbw      = array( NA, dim = dim.envr)
      ts.leaf.gsw      = array( NA, dim = dim.envr)
      ts.cba.light     = array( NA, dim = dim.envr)
      ts.cba.moist     = array( NA, dim = dim.envr)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop over all times, and retrieve the data.                                    #
      #------------------------------------------------------------------------------------#
      for (u in 1:nfullcyc){
         now.month = nummonths(sta$when[1])
         now.year  = i1stfull - dyear[2] + (u-1) * ncyc
         year.use  = act.census.yeara - dyear[2]
         for (m in 1:n.months){
            now.month   = (now.month %% 12) + 1
            now.year    = now.year + as.integer(now.month == 1)
            year.use    = year.use + as.integer(now.month == 1)
            cat("Cycle",u,"; ED-2.1: ",paste(mon2mmm(now.month,cap1=T),now.year,sep="-")
                         ,"; Census: ",paste(act.census.year[census.idx[m]]),"\n")

            #----- Build the file name. ---------------------------------------------------#
            cmonth = sprintf("%2.2i",now.month)
            cyear  = sprintf("%2.2i",now.year )
            myfile = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5",sep="")
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    muse is the m that should be used for agb, basal area, and LAI.           #
            #------------------------------------------------------------------------------#
            if (m == 1){
               muse = 1
            }else{
               if (census.idx[m] != census.idx[m-1]){
                  muse = m
               }#end if
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Read data if the file exists.                                            #
            #------------------------------------------------------------------------------#
            if (file.exists(myfile)){

               #----- Read data and close connection immediately after. -------------------#
               #cat("     * Reading ",basename(myfile),"...","\n")
               mymont = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
               #---------------------------------------------------------------------------#


               #---- Read soil information if it hasn't been read yet. --------------------#
               if (read.soil){
                  read.soil  = FALSE
                  nzg        = mymont$NZG
                  isoilflg   = mymont$ISOILFLG
                  slz        = mymont$SLZ
                  slxsand    = mymont$SLXSAND
                  slxclay    = mymont$SLXCLAY
                  ntext      = mymont$NTEXT.SOIL[nzg]
                  soil       = soil.params(ntext,isoilflg,slxsand,slxclay)
                  dslz       = diff(c(slz,0))
                  soil.depth = rev(cumsum(rev(dslz)))
                  soil.dry   = rev(cumsum(rev(soil$soilcp * wdns * dslz)))
                  soil.poro  = rev(cumsum(rev(soil$slmsts * wdns * dslz)))


                  #----- Find the layers we care about. -----------------------------------#
                  sel        = slz < slz.min
                  if (any(sel)){
                     ka      = which.max(slz[sel])
                  }else{
                     ka      = 1
                  }#end if
                  kz         = nzg
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Save the year and the census index.                                  #
               #---------------------------------------------------------------------------#
               ts.census.year[m,u] = year.use
               ts.census.idx [m,u] = census.idx [m]
               #---------------------------------------------------------------------------#


               #---- Read in the site-level area. -----------------------------------------#
               areasi     = mymont$AREA.SI
               npatches   = mymont$SIPA.N
               #---------------------------------------------------------------------------#


               #----- Read a few patch-level variables. -----------------------------------#
               areapa     = mymont$AREA * rep(areasi,times=npatches)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Read in the soil moisture, and find the equivalent matric potential.  #
               #---------------------------------------------------------------------------#
               soil.water       = rev(cumsum(rev(mymont$MMEAN.SOIL.WATER.PY * wdns * dslz)))
               soil.moist.avg   = soil.water / (wdns * soil.depth)
               ts.paw     [m,u] = ( ( soil.water[ka] - soil.dry [ka] )
                                  / ( soil.poro [ka] - soil.dry [ka] ) )
               ts.rshort  [m,u] = mymont$MMEAN.ATM.RSHORT.PY
               ts.fs.open [m,u] = mymont$MMEAN.FS.OPEN.PY
               ts.smpot   [m,u] = ( - smoist2mpot(smoist=soil.moist.avg[ka],mysoil=soil)
                                  * 0.001 * grav )
               ts.atm.vpd [m,u] = mymont$MMEAN.ATM.VPDEF.PY
               ts.leaf.vpd[m,u] = mymont$MMEAN.LEAF.VPDEF.PY
               ts.leaf.gbw[m,u] = mymont$MMEAN.LEAF.GBW.PY
               ts.leaf.gsw[m,u] = mymont$MMEAN.LEAF.GSW.PY
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Read the cohort-level variables.  Because empty patchs do exist       #
               # (deserts), we must check whether there is any cohort to be read.  If not, #
               # assign NA to all variables.                                               #
               #---------------------------------------------------------------------------#
               ncohorts   = mymont$PACO.N
               if (any (ncohorts > 0)){
                  #----- Make a cohort-level area. ----------------------------------------#
                  areaconow    = rep(areapa,times=ncohorts)
                  #------------------------------------------------------------------------#


                  #----- Define the DBH classes. ------------------------------------------#
                  dbhconow     = mymont$DBH
                  dbhcut       = cut(dbhconow,breaks=sta$dbh.breaks)
                  dbhlevs      = levels(dbhcut)
                  dbhfac       = match(dbhcut,dbhlevs)
                  n.dbh        = length(dbhlevs)
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Load the other cohort-level variables of interest.                 #
                  #------------------------------------------------------------------------#
                  pftconow        = mymont$PFT
                  nplantconow     = mymont$NPLANT
                  agbconow        = mymont$AGB.CO
                  baconow         = mymont$BA.CO
                  laiconow        = mymont$LAI.CO
                  mortconow       = rowSums(mymont$MMEAN.MORT.RATE.CO)
                  ddmortconow     = mymont$MMEAN.MORT.RATE.CO[,2]
                  dimortconow     = mortconow - ddmortconow
                  recruitconow    = mymont$RECRUIT.DBH
                  censtatusconow  = mymont$CENSUS.STATUS
                  dlndbhdtconow   = mymont$DLNDBH.DT
                  dlnagbdtconow   = mymont$DLNAGB.DT
                  dlnbadtconow    = mymont$DLNBA.DT
                  cbalightconow   = rowMeans(mymont$CB.LIGHTMAX[,1:12])
                  cbamoistconow   = rowMeans(mymont$CB.MOISTMAX[,1:12])
                  #------------------------------------------------------------------------#
               }else{
                  areaconow       = NA
                  dbhconow        = NA
                  pftconow        = NA
                  nplantconow     = NA
                  agbconow        = NA
                  baconow         = NA
                  laiconow        = NA
                  mortconow       = NA
                  ddmortconow     = NA
                  dimortconow     = NA
                  recruitconow    = NA
                  censtatusconow  = NA
                  dlndbhdtconow   = NA
                  dlnagbdtconow   = NA
                  dlnbadtconow    = NA
                  cbalightconow   = NA
                  cbamoistconow   = NA
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     The following variables are used to scale "intensive" properties      #
               # (whatever/plant) to "extensive" (whatever/m2).  Sometimes it may be used  #
               # to build weighted averages.                                               #
               #---------------------------------------------------------------------------#
               w.nplant = nplantconow * areaconow
               w.lai    = laiconow    * areaconow
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Find the polygon-level carbon balance.                               #
               #---------------------------------------------------------------------------#
               if (all(!is.na(nplantconow))){
                  ts.cba.light[m,u] = sum( w.nplant * cbalightconow )
                  ts.cba.moist[m,u] = sum( w.nplant * cbamoistconow )
               }else{
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over all PFTs.                                                   #
               #---------------------------------------------------------------------------#
               for (p in 1:(npft+1)){
                  #------------------------------------------------------------------------#
                  #     Select the subsample.                                              #
                  #------------------------------------------------------------------------#
                  if (all(is.na(pftconow))){
                     s.dbh  = rep(FALSE,times=length(pftconow))
                     s.cs2  = rep(FALSE,times=length(pftconow))
                  }else if (p <= npft){
                     s.dbh  = pftconow == p & censtatusconow >  0
                     s.cs2  = pftconow == p & censtatusconow == 2
                  }else{
                     s.dbh  = censtatusconow >  0
                     s.cs2  = censtatusconow == 2
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #    Find the AGB, LAI, and BA if this is the first time, otherwise,     #
                  # copy from the first time.  This is because an actual census would not  #
                  # have the information.                                                  #
                  #------------------------------------------------------------------------#
                  if (m == muse && any(s.dbh)){
                     ts.n.plot   [p,m,u] = sum( w.nplant[s.dbh]                   )
                     ts.agb.plot [p,m,u] = sum( w.nplant[s.dbh] * agbconow[s.dbh] )
                     ts.ba.plot  [p,m,u] = sum( w.nplant[s.dbh] * baconow [s.dbh] )
                     ts.lai.plot [p,m,u] = sum( w.lai   [s.dbh]                   )
                  }else{
                     ts.n.plot   [p,m,u] = ts.n.plot   [p,muse,u]
                     ts.agb.plot [p,m,u] = ts.agb.plot [p,muse,u]
                     ts.ba.plot  [p,m,u] = ts.ba.plot  [p,muse,u]
                     ts.lai.plot [p,m,u] = ts.lai.plot [p,muse,u]
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Find the PFT-level mortality rates.                                #
                  #------------------------------------------------------------------------#
                  if (any(s.cs2)){

                     #---- This is the number of survivors. -------------------------------#
                     survivor                = sum( w.nplant       [s.cs2] )
                     previous                = sum( w.nplant       [s.cs2]
                                                  * exp(mortconow  [s.cs2]))
                     ts.mort.plot$n  [p,m,u] = log( previous / survivor )

                     survivor                = sum( w.nplant       [s.cs2])
                     previous                = sum( w.nplant       [s.cs2]
                                                  * exp(dimortconow[s.cs2]))
                     ts.dimort.plot$n[p,m,u] = log( previous / survivor )

                     survivor                = sum( w.nplant       [s.cs2])
                     previous                = sum( w.nplant       [s.cs2]
                                                  * exp(ddmortconow[s.cs2]))
                     ts.ddmort.plot$n[p,m,u] = log( previous / survivor )
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Find the PFT-level recruitment rates.                              #
                  #------------------------------------------------------------------------#
                  if (any(s.dbh) & any(s.cs2)){
                     #---- This is the number of survivors. -------------------------------#
                     population              = sum(w.nplant[s.dbh])
                     established             = sum(w.nplant[s.cs2])
                     ts.recr.plot$n  [p,m,u] = log( population / established) / 12.0
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Growth rates are found only for established cohorts.               #
                  #------------------------------------------------------------------------#
                  if (any(s.cs2)){
                     #----- Estimate the median of this population. -----------------------#
                     dlndbhdt        = dlndbhdtconow [s.cs2]
                     dlnagbdt        = dlnagbdtconow [s.cs2]
                     dlnbadt         = dlnbadtconow  [s.cs2]
                     wgt             = w.nplant      [s.cs2]
                     median.dlndbhdt = weighted.quantile( x  = dlndbhdtconow [s.cs2]
                                                        , w  = w.nplant      [s.cs2]
                                                        , qu = 0.50
                                                        )#end 
                     median.dlnagbdt = weighted.quantile( x  = dlnagbdtconow [s.cs2]
                                                        , w  = w.nplant      [s.cs2]
                                                        , qu = 0.50
                                                        )#end 
                     median.dlnbadt  = weighted.quantile( x  = dlnbadtconow  [s.cs2]
                                                        , w  = w.nplant      [s.cs2]
                                                        , qu = 0.50
                                                        )#end 
                     #---------------------------------------------------------------------#

                     #---------------------------------------------------------------------#
                     #     Copy the results to the arrays.                                 #
                     #---------------------------------------------------------------------#
                     ts.growth.plot$dbh[p,m,u] = median.dlndbhdt$q
                     ts.growth.plot$agb[p,m,u] = median.dlnagbdt$q
                     ts.growth.plot$ba [p,m,u] = median.dlnbadt$q
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Build the size (DBH) structure arrays.                                #
               #---------------------------------------------------------------------------#
               for (d in 1:n.dbh){
                  if (all(is.na(dbhfac))){
                     this.dbh  = rep(FALSE,times=length(dbhfac))
                  }else{
                     this.dbh  = dbhfac == d
                  }#end if
                  for (p in 1:(npft+1)){
                     #---------------------------------------------------------------------#
                     #      Select the cohorts that will be used here.                     #
                     #---------------------------------------------------------------------#
                     if (p <= npft){
                        s.pft  = pftconow == p
                        s.cs2  = s.pft & this.dbh & censtatusconow == 2
                        s.dbh  = s.pft & this.dbh & censtatusconow  > 0
                     }else{
                        s.cs2  = this.dbh & censtatusconow == 2
                        s.dbh  = this.dbh & censtatusconow  > 0
                     }#end if
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #    Find the AGB, LAI, and BA if this is the first time, otherwise,  #
                     # copy from the first time.  This is because an actual census would   #
                     # not have the information.                                           #
                     #---------------------------------------------------------------------#
                     if (m == muse && any(s.dbh)){
                        ts.n.size   [p,d,m,u] = sum( w.nplant[s.dbh]                   )
                        ts.agb.size [p,d,m,u] = sum( w.nplant[s.dbh] * agbconow[s.dbh] )
                        ts.ba.size  [p,d,m,u] = sum( w.nplant[s.dbh] * baconow [s.dbh] )
                        ts.lai.size [p,d,m,u] = sum( w.lai   [s.dbh]                   )
                     }else{
                        ts.n.size   [p,d,m,u] = ts.n.size   [p,d,muse,u]
                        ts.agb.size [p,d,m,u] = ts.agb.size [p,d,muse,u]
                        ts.ba.size  [p,d,m,u] = ts.ba.size  [p,d,muse,u]
                        ts.lai.size [p,d,m,u] = ts.lai.size [p,d,muse,u]
                     }#end if
                     #---------------------------------------------------------------------#

                     if (any(s.cs2)){
                        #---- This is the number of survivors and living before. ----------#
                        survivor                   = sum( w.nplant        [s.cs2] )
                        previous                   = sum( w.nplant        [s.cs2]
                                                        * exp(mortconow   [s.cs2]))
                        ts.mort.size$n   [p,d,m,u] = log( previous / survivor )

                        survivor                   = sum( w.nplant        [s.cs2])
                        previous                   = sum( w.nplant        [s.cs2] 
                                                        * exp(dimortconow [s.cs2]))
                        ts.dimort.size$n [p,d,m,u] = log( previous / survivor )

                        survivor                   = sum( w.nplant        [s.cs2])
                        previous                   = sum( w.nplant        [s.cs2] 
                                                        * exp(ddmortconow [s.cs2]))
                        ts.ddmort.size$n [p,d,m,u] = log( previous / survivor )
                     }#end if
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Growth rates are found only for established cohorts.            #
                     #---------------------------------------------------------------------#
                     if (any(s.cs2)){
                        #----- Estimate the median of this population. --------------------#
                        dlndbhdt        = dlndbhdtconow [s.cs2]
                        dlnagbdt        = dlnagbdtconow [s.cs2]
                        dlnbadt         = dlnbadtconow  [s.cs2]
                        wgt             = w.nplant      [s.cs2]
                        median.dlndbhdt = weighted.quantile( x  = dlndbhdtconow [s.cs2]
                                                           , w  = w.nplant      [s.cs2]
                                                           , qu = 0.50
                                                           )#end 
                        median.dlnagbdt = weighted.quantile( x  = dlnagbdtconow [s.cs2]
                                                           , w  = w.nplant      [s.cs2]
                                                           , qu = 0.50
                                                           )#end 
                        median.dlnbadt  = weighted.quantile( x  = dlnbadtconow  [s.cs2]
                                                           , w  = w.nplant      [s.cs2]
                                                           , qu = 0.50
                                                           )#end 
                        #------------------------------------------------------------------#

                        #------------------------------------------------------------------#
                        #     Copy the results to the arrays.                              #
                        #------------------------------------------------------------------#
                        ts.growth.size$dbh[p,d,m,u] = median.dlndbhdt$q
                        ts.growth.size$agb[p,d,m,u] = median.dlnagbdt$q
                        ts.growth.size$ba [p,d,m,u] = median.dlnbadt$q
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for PFT
                  #------------------------------------------------------------------------#
               }#end for DBH
               #---------------------------------------------------------------------------#
            }else{
               cat("     * ",basename(myfile)," wasn't found, skipping it...","\n")
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (m in 1:n.months)
         #---------------------------------------------------------------------------------#
      }#end for (u in 1:nfullcyc)
      #====================================================================================#
      #====================================================================================#




      #====================================================================================#
      #====================================================================================#
      #     Make the parameter space.                                                      #
      #------------------------------------------------------------------------------------#
      cat("   - Creating the parameter space...","\n")
         #----- Make the point and legend name sequence. ----------------------------------#
         yeara            = c(qapply(X=ts.census.year,INDEX=census.idx,DIM=1
                                    ,FUN=min,na.rm=TRUE))
         yearz            = c(qapply(X=ts.census.year,INDEX=census.idx,DIM=1
                                    ,FUN=max,na.rm=TRUE))
         pspace$pch       = eft.pch[match(yearz,eft.year)]
         pspace$leg.label = paste(sort(unique(yeara[is.finite(yeara)]))
                                 ,sort(unique(yearz[is.finite(yearz)])),sep="-")
         pspace$leg.pch   = eft.pch[match(sort(unique(yearz[is.finite(yearz)])),eft.year)]
         #---------------------------------------------------------------------------------#



         #----- Environment variables. ----------------------------------------------------#
         pspace$rshort           = c(qapply( X     = ts.rshort
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         pspace$fs.open          = c(qapply( X     = ts.fs.open
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.20
                                           , na.rm = TRUE
                                           , lower = TRUE
                                           ))
         pspace$paw              = c(qapply( X     = ts.paw
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.20
                                           , na.rm = TRUE
                                           , lower = TRUE
                                           ))
         pspace$smpot            = c(qapply( X     = ts.smpot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         pspace$atm.vpd          = c(qapply( X     = ts.atm.vpd
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         pspace$leaf.vpd         = c(qapply( X     = ts.leaf.vpd
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         pspace$leaf.gsw        = c(qapply( X     = ts.leaf.gsw
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         pspace$leaf.gbw        = c(qapply( X     = ts.leaf.gbw
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         pspace$cba.light       = c(qapply( X     = ts.cba.light
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         pspace$cba.moist       = c(qapply( X     = ts.cba.moist
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = qu.mean
                                           , p     = 0.80
                                           , na.rm = TRUE
                                           , lower = FALSE
                                           ))
         #---------------------------------------------------------------------------------#



         #----- Plot-level variables. -----------------------------------------------------#
         pspace$lai.plot.mean     = c(qapply( X     = ts.lai.plot       [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean 
                                            , na.rm = TRUE
                                            ))
         pspace$ba.plot.mean      = c(qapply( X     = ts.ba.plot        [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean 
                                            , na.rm = TRUE
                                            ))
         pspace$agb.plot.mean     = c(qapply( X     = ts.agb.plot       [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean 
                                            , na.rm = TRUE
                                            ))
         pspace$recr.mean         = c(qapply( X     = ts.recr.plot$n    [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean 
                                            , na.rm = TRUE
                                            ))
         pspace$mort.plot.mean    = c(qapply( X     = ts.mort.plot$n    [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ))
         pspace$dimort.plot.mean  = c(qapply( X     = ts.dimort.plot$n  [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ))
         pspace$ddmort.plot.mean  = c(qapply( X     = ts.ddmort.plot$n  [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ))
         pspace$growdbh.plot.mean = c(qapply( X     = ts.growth.plot$dbh[npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ))
         pspace$growagb.plot.mean = c(qapply( X     = ts.growth.plot$agb[npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ))
         pspace$growba.plot.mean  = c(qapply( X     = ts.growth.plot$ba [npft+1,,]
                                            , INDEX = census.idx
                                            , DIM   = 1
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ))
         #---------------------------------------------------------------------------------#



         #----- Size-level variables. -----------------------------------------------------#
         pspace$lai.size.mean     = matrix(qapply( X     = ts.lai.size   [npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$ba.size.mean      = matrix(qapply( X     = ts.ba.size    [npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$agb.size.mean     = matrix(qapply( X     = ts.agb.size   [npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$mort.size.mean    = matrix(qapply( X     = ts.mort.size$n[npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$dimort.size.mean  = matrix(qapply( X     = ts.dimort.size$n  [npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$ddmort.size.mean  = matrix(qapply( X     = ts.ddmort.size$n  [npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$growdbh.size.mean = matrix(qapply( X     = ts.growth.size$dbh[npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$growagb.size.mean = matrix(qapply( X     = ts.growth.size$agb[npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         pspace$growba.size.mean  = matrix(qapply( X     = ts.growth.size$ba [npft+1,,,]
                                                 , INDEX = census.idx
                                                 , DIM   = 2
                                                 , FUN   = mean
                                                 , na.rm = TRUE
                                                 )
                                          , nrow = n.dbh
                                          , ncol = nfullcyc*(n.census-1)
                                          )#end matrix
         #---------------------------------------------------------------------------------#
      #====================================================================================#
      #====================================================================================#






      #====================================================================================#
      #====================================================================================#
      #------------------------------------------------------------------------------------#
      #     Remove the 4th. dimension for the means.                                       #
      #------------------------------------------------------------------------------------#
      cat("   - Averaging the census intervals...","\n")


      #----- This function will be used by the bootstrap. ---------------------------------#
      mean.fun = function(x,idx) mean(x[idx],na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot-level recruitment.                                                        #
      #------------------------------------------------------------------------------------#
      ms.recr.plot   = list()
      ms.mort.plot   = list()
      ms.ddmort.plot = list()
      ms.dimort.plot = list()
      ms.growth.plot = list()
      ms.mort.size   = list()
      ms.ddmort.size = list()
      ms.dimort.size = list()
      ms.growth.size = list()
      #------------------------------------------------------------------------------------#



      #----- Plot-level recruitment. ------------------------------------------------------#
      cat("     * Recruitment...","\n")
      for (v in 1:nrecr.vars){
         v.now  = recr.vars[v]
         ms.now = list()

         recr.rates = c("recr")
         n.rates    = length(recr.rates)

         #---------------------------------------------------------------------------------#
         #     Loop over all rates.                                                        #
         #---------------------------------------------------------------------------------#
         for (r in 1:n.rates){
            ts.this.plot = get(paste("ts",recr.rates[r],"plot",sep="."))
            ms.this.plot = get(paste("ms",recr.rates[r],"plot",sep="."))

            #------------------------------------------------------------------------------#
            #     Find the mean for all PFTs.                                              #
            #------------------------------------------------------------------------------#
            ms.mean.plot = array(NA,dim=c(npft+1,n.census))
            ms.q025.plot = array(NA,dim=c(npft+1,n.census))
            ms.q975.plot = array(NA,dim=c(npft+1,n.census))

            for (p in 1:(npft+1)){
               for (i in 2:n.census){
                  i.sel = census.idx == i
                  ts.plot.now       = c(ts.this.plot[[v.now]][p,i.sel,])
                  ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                  if (any(is.finite(ts.plot.now))){
                     boot.now = boot   (data=ts.plot.now,statistic=mean.fun,R=n.boot)
                     ci.now   = boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                     if (length(ci.now$percent) == 5){
                        ms.q025.plot[p,i] = ci.now$percent[4]
                        ms.q975.plot[p,i] = ci.now$percent[5]
                     }else{
                        warning("Failed using bootstrap...")
                     }#end if
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#



            #----- Convert rates fraction rates. ------------------------------------------#
            ms.mean.plot = exp(ms.mean.plot) - 1.
            ms.q025.plot = exp(ms.q025.plot) - 1.
            ms.q975.plot = exp(ms.q975.plot) - 1.
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Save plot and size.                                                      #
            #------------------------------------------------------------------------------#
            ms.this.plot[[v.now]] = list( mean = ms.mean.plot
                                        , q025 = ms.q025.plot
                                        , q975 = ms.q975.plot
                                        )#end list
            #------------------------------------------------------------------------------#
            dummy = assign(paste("ms",recr.rates[r],"plot",sep="."), ms.this.plot)
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#


      #----- Mortality rates. -------------------------------------------------------------#
      cat("     * Mortality...","\n")
      for (v in 1:nmort.vars){
         v.now  = mort.vars[v]
         ms.now = list()

         mort.rates = c("mort","ddmort","dimort")
         n.rates    = length(mort.rates)

         #---------------------------------------------------------------------------------#
         #     Loop over all rates.                                                        #
         #---------------------------------------------------------------------------------#
         for (r in 1:n.rates){
            ts.this.plot = get(paste("ts",mort.rates[r],"plot",sep="."))
            ts.this.size = get(paste("ts",mort.rates[r],"size",sep="."))
            ms.this.plot = get(paste("ms",mort.rates[r],"plot",sep="."))
            ms.this.size = get(paste("ms",mort.rates[r],"size",sep="."))


            #------------------------------------------------------------------------------#
            #     Find the mean for all PFTs.                                              #
            #------------------------------------------------------------------------------#
            ms.mean.plot = array(NA,dim=c(npft+1,n.census))
            ms.q025.plot = array(NA,dim=c(npft+1,n.census))
            ms.q975.plot = array(NA,dim=c(npft+1,n.census))
            ms.mean.size = array(NA,dim=c(npft+1,n.dbh,n.census))
            ms.q025.size = array(NA,dim=c(npft+1,n.dbh,n.census))
            ms.q975.size = array(NA,dim=c(npft+1,n.dbh,n.census))

            for (p in 1:(npft+1)){
               for (i in 2:n.census){
                  i.sel = census.idx == i
                  ts.plot.now       = c(ts.this.plot[[v.now]][p,i.sel,])
                  ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                  if (any(is.finite(ts.plot.now))){
                     boot.now = boot   (data=ts.plot.now,statistic=mean.fun,R=n.boot)
                     ci.now   = boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                     if (length(ci.now$percent) == 5){
                        ms.q025.plot[p,i] = ci.now$percent[4]
                        ms.q975.plot[p,i] = ci.now$percent[5]
                     }else{
                        warning("Failed using bootstrap...")
                     }#end if
                  }#end if
                  #------------------------------------------------------------------------#

                  #------------------------------------------------------------------------#
                  for (d in 1:n.dbh){
                     ts.size.now         = c(ts.this.size[[v.now]][p,d,i.sel,])
                     ms.mean.size[p,d,i] = mean(ts.size.now,na.rm=TRUE) 
                     if (any(is.finite(ts.size.now))){
                        boot.now = boot   (data=ts.size.now,statistic=mean.fun,R=n.boot)
                        ci.now   = boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                        if (length(ci.now$percent) == 5){
                           ms.q025.size[p,d,i] = ci.now$percent[4]
                           ms.q975.size[p,d,i] = ci.now$percent[5]
                        }else{
                           warning("Failed using bootstrap...")
                        }#end if
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#



            #----- Convert rates fraction rates. ------------------------------------------#
            ms.mean.plot = 1. - exp( - ms.mean.plot)
            ms.q025.plot = 1. - exp( - ms.q025.plot)
            ms.q975.plot = 1. - exp( - ms.q975.plot)
            ms.mean.size = 1. - exp( - ms.mean.size)
            ms.q025.size = 1. - exp( - ms.q025.size)
            ms.q975.size = 1. - exp( - ms.q975.size)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Save plot and size.                                                      #
            #------------------------------------------------------------------------------#
            ms.this.plot[[v.now]] = list( mean = ms.mean.plot
                                        , q025 = ms.q025.plot
                                        , q975 = ms.q975.plot
                                        )#end list
            ms.this.size[[v.now]] = list( mean = ms.mean.size
                                        , q025 = ms.q025.size
                                        , q975 = ms.q975.size
                                        )#end list
            #------------------------------------------------------------------------------#
            dummy = assign(paste("ms",mort.rates[r],"plot",sep="."), ms.this.plot)
            dummy = assign(paste("ms",mort.rates[r],"size",sep="."), ms.this.size)
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#

      
      #----- Growth rates. ----------------------------------------------------------------#
      cat("     * Growth...","\n")
      for (v in 1:ngrowth.vars){
         v.now  = growth.vars[v]
         ms.now = list()

         growth.rates = c("growth")
         n.rates    = length(growth.rates)

         #---------------------------------------------------------------------------------#
         #     Loop over all rates.                                                        #
         #---------------------------------------------------------------------------------#
         for (r in 1:n.rates){
            ts.this.plot = get(paste("ts",growth.rates[r],"plot",sep="."))
            ts.this.size = get(paste("ts",growth.rates[r],"size",sep="."))
            ms.this.plot = get(paste("ms",growth.rates[r],"plot",sep="."))
            ms.this.size = get(paste("ms",growth.rates[r],"size",sep="."))

            #------------------------------------------------------------------------------#
            #     Find the mean for all PFTs.                                              #
            #------------------------------------------------------------------------------#
            ms.mean.plot = array(NA,dim=c(npft+1,n.census))
            ms.q025.plot = array(NA,dim=c(npft+1,n.census))
            ms.q975.plot = array(NA,dim=c(npft+1,n.census))
            ms.mean.size = array(NA,dim=c(npft+1,n.dbh,n.census))
            ms.q025.size = array(NA,dim=c(npft+1,n.dbh,n.census))
            ms.q975.size = array(NA,dim=c(npft+1,n.dbh,n.census))

            for (p in 1:(npft+1)){
               for (i in 2:n.census){
                  i.sel = census.idx == i
                  ts.plot.now       = c(ts.this.plot[[v.now]][p,i.sel,])
                  ms.mean.plot[p,i] = mean(ts.plot.now,na.rm=TRUE) 
                  if (any(is.finite(ts.plot.now))){
                     boot.now = boot   (data=ts.plot.now,statistic=mean.fun,R=n.boot)
                     ci.now   = boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                     if (length(ci.now$percent) == 5){
                        ms.q025.plot[p,i] = ci.now$percent[4]
                        ms.q975.plot[p,i] = ci.now$percent[5]
                     }else{
                        warning("Failed using bootstrap...")
                     }#end if
                  }#end if
                  #------------------------------------------------------------------------#

                  #------------------------------------------------------------------------#
                  for (d in 1:n.dbh){
                     ts.size.now         = c(ts.this.size[[v.now]][p,d,i.sel,])
                     ms.mean.size[p,d,i] = mean(ts.size.now,na.rm=TRUE) 
                     if (any(is.finite(ts.size.now))){
                        boot.now = boot   (data=ts.size.now,statistic=mean.fun,R=n.boot)
                        ci.now   = boot.ci(boot.out=boot.now,conf=0.95,type="perc")
                        if (length(ci.now$percent) == 5){
                           ms.q025.size[p,d,i] = ci.now$percent[4]
                           ms.q975.size[p,d,i] = ci.now$percent[5]
                        }else{
                           warning("Failed using bootstrap...")
                        }#end if
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Save plot and size.                                                      #
            #------------------------------------------------------------------------------#
            ms.this.plot[[v.now]] = list( mean = ms.mean.plot
                                        , q025 = ms.q025.plot
                                        , q975 = ms.q975.plot
                                        )#end list
            ms.this.size[[v.now]] = list( mean = ms.mean.size
                                        , q025 = ms.q025.size
                                        , q975 = ms.q975.size
                                        )#end list
            #------------------------------------------------------------------------------#
            dummy = assign(paste("ms",growth.rates[r],"plot",sep="."), ms.this.plot)
            dummy = assign(paste("ms",growth.rates[r],"size",sep="."), ms.this.size)
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the average rates for the census period using log-normal.                 #
      #------------------------------------------------------------------------------------#
      mypfts = sort(match(unique(names(sta$classes)),pft$name))
      npfts  = length(mypfts)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Retrieve the factor, classes and wood density used for the observations.  We  #
      # can't switch the factors between observation and statistics because we use         #
      # observations to drive the statistics.                                              #
      #------------------------------------------------------------------------------------#
      cat("   - Finding the average rates...","\n")
      ed2           = list()
      ed2$when      = sta$when
      ed2$taxon     = sta$taxon
      ed2$classes   = sta$classes
      nfac          = length(ed2$classes)
      ed2$wood.dens = sta$wood.dens
      ed2$dtime     = sta$dtime
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Recruitment rates.  Only global values can be found.                          #
      #------------------------------------------------------------------------------------#
      for (r in 1:npratets){
         #---------------------------------------------------------------------------------#
         #     Load the rate information.                                                  #
         #---------------------------------------------------------------------------------#
         this.rate = pratets[[ r]]
         ed2.rate  = this.rate$ed2.rate
         sta.rate  = this.rate$sta.rate
         sizetoo   = this.rate$sizetoo
         desc.rate = this.rate$desc.rate
         indiv     = this.rate$indiv
         cat(" - Compounding the ",desc.rate," tables...","\n")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Load plot-level and size-level data.                                       #
         #---------------------------------------------------------------------------------#
         ed2.plot   = paste(ed2.rate,"plot",sep=".")
         sta.plot   = paste(sta.rate,"plot",sep=".")
         ed2.size   = paste(ed2.rate,"size",sep=".")
         sta.size   = paste(sta.rate,"size",sep=".")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Find how many individuals to retrieve.                                    #
         #---------------------------------------------------------------------------------#
         nindiv = length(indiv)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop over the different types of individuals.                               #
         #---------------------------------------------------------------------------------#
         ed2[[ed2.plot]] = list()
         for (v in 1:nindiv){
            #----- Set up the individuals. ------------------------------------------------#
            v.now   = indiv[v]
            yyy     = 2:n.census
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot-level with all cycles.                                              #
            #------------------------------------------------------------------------------#
            ms.plot   = get(paste("ms",ed2.rate,"plot",sep="."))[[v.now]]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Create the plot-level structure.  Here the median and quantiles have a   #
            # somewhat different number, based on the number of cycles we run (median      #
            # amongst each cycle).                                                         #
            #------------------------------------------------------------------------------#
            ed2[[ed2.plot]][[v.now]]               = list()
            #----- "Borrow" the structure from the sta counterpart. ------------------------#
            ed2[[ed2.plot]][[v.now]]$global        = NA * sta[[sta.plot]][[v.now]]$global
            ed2[[ed2.plot]][[v.now]]$median        = NA * sta[[sta.plot]][[v.now]]$median
            ed2[[ed2.plot]][[v.now]]$q025          = NA * sta[[sta.plot]][[v.now]]$q025
            ed2[[ed2.plot]][[v.now]]$q975          = NA * sta[[sta.plot]][[v.now]]$q975
            #----- Save the global variables. ---------------------------------------------#
            ed2[[ed2.plot]][[v.now]]$global[4,yyy] = ms.plot$mean[npft+1,yyy]
            ed2[[ed2.plot]][[v.now]]$global[5,yyy] = ms.plot$q025[npft+1,yyy]
            ed2[[ed2.plot]][[v.now]]$global[6,yyy] = ms.plot$q975[npft+1,yyy]
            #----- Save the PFT statistics. -----------------------------------------------#
            ed2[[ed2.plot]][[v.now]]$median[ ,yyy] = ms.plot$mean[mypfts,yyy]
            ed2[[ed2.plot]][[v.now]]$q025  [ ,yyy] = ms.plot$q025[mypfts,yyy]
            ed2[[ed2.plot]][[v.now]]$q975  [ ,yyy] = ms.plot$q975[mypfts,yyy]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Now the size-dependent variables, if this is a size-dependent rate.       #
            #------------------------------------------------------------------------------#
            if (sizetoo){

               #---------------------------------------------------------------------------#
               #     Plot-level with all cycles.                                           #
               #---------------------------------------------------------------------------#
               ms.size   = get(paste("ms",ed2.rate,"size",sep="."))[[v.now]]
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Create the plot-level structure.  Here the median and quantiles have  #
               # a somewhat different number, based on the number of cycles we run (median #
               # amongst each cycle).                                                      #
               #---------------------------------------------------------------------------#
               ed2[[ed2.size]][[v.now]]        = list()
               #----- "Borrow" the structure from the sta counterpart. --------------------#
               ed2[[ed2.size]][[v.now]]$global = NA * sta[[sta.size]][[v.now]]$global
               ed2[[ed2.size]][[v.now]]$median = NA * sta[[sta.size]][[v.now]]$median
               ed2[[ed2.size]][[v.now]]$q025   = NA * sta[[sta.size]][[v.now]]$q025
               ed2[[ed2.size]][[v.now]]$q975   = NA * sta[[sta.size]][[v.now]]$q975
               #----- Save the global variables. ------------------------------------------#
               ed2[[ed2.size]][[v.now]]$global[4,,yyy] = ms.size$mean[npft+1,,yyy]
               ed2[[ed2.size]][[v.now]]$global[5,,yyy] = ms.size$q025[npft+1,,yyy]
               ed2[[ed2.size]][[v.now]]$global[6,,yyy] = ms.size$q975[npft+1,,yyy]
               #----- Save the PFT statistics. --------------------------------------------#
               ed2[[ed2.size]][[v.now]]$median[ ,,yyy] = ms.size$mean[mypfts,,yyy]
               ed2[[ed2.size]][[v.now]]$q025  [ ,,yyy] = ms.size$q025[mypfts,,yyy]
               ed2[[ed2.size]][[v.now]]$q975  [ ,,yyy] = ms.size$q975[mypfts,,yyy]
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Make the directories.                                                          #
      #------------------------------------------------------------------------------------#
      outplot   = paste(outpref,"census_plot"  ,sep="/")
      outsize   = paste(outpref,"census_size"  ,sep="/")
      if (! file.exists(outplot)) dir.create(outplot)
      if (! file.exists(outsize)) dir.create(outsize)
      #------------------------------------------------------------------------------------#









      #====================================================================================#
      #====================================================================================#
      #  6.  Plot rates as function of size.                                               #
      #====================================================================================#
      #====================================================================================#
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



         #----- Create a directory for this type of plot. ---------------------------------#
         outrate   = paste(outsize,ed2.rate,sep="/")
         if (! file.exists(outrate)) dir.create(outrate)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Loop over all possible types of population count.                            #
         #---------------------------------------------------------------------------------#
         for (i in 1:nindiv){
            cat("  - Size level: ",desc.indiv[i],"...","\n")


            #----- Build the rate name. ---------------------------------------------------#
            cat(" + Plotting size-dependent ",desc.rate,"...","\n")
            ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
            sta.rate = paste(this.plot$sta.rate,"size",sep=".")
            #------------------------------------------------------------------------------#



            #----- Create path for this individual. ---------------------------------------#
            outindiv   = paste(outrate,indiv[i],sep="/")
            if (! file.exists(outindiv)) dir.create(outindiv)
            #------------------------------------------------------------------------------#


            #----- Load the modelled rates. -----------------------------------------------#
            sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
            sta.median   = 100. * sta.mod[4,,]
            sta.q025     = 100. * sta.mod[5,,]
            sta.q975     = 100. * sta.mod[6,,]
            ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
            ed2.median   = 100. * ed2.mod[4,,]
            ed2.q025     = 100. * ed2.mod[5,,]
            ed2.q975     = 100. * ed2.mod[6,,]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Define a nice configuration for the multiple panels.                    #
            #------------------------------------------------------------------------------#
            lo.box = pretty.box(n=n.census-1,horizontal=TRUE)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all formats.                                                  #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open the file or the plot window. -----------------------------------#
               fichier = paste(outindiv,"/yrsize-",ed2.rate,"-",indiv[i],".",outform[o]
                              ,sep="")
               if(outform[o] == "x11"){
                  X11(width=wide.size$width,height=wide.size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=wide.size$width*depth
                     ,height=wide.size$height*depth,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=wide.size$width,height=wide.size$height
                            ,pointsize=ptsz,paper=wide.size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=wide.size$width,height=wide.size$height,pointsize=ptsz
                     ,paper=wide.size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Split the window into several smaller windows.                        #
               #---------------------------------------------------------------------------#
               par(par.user)
               par.orig = par(no.readonly = TRUE)
               par(oma = c(0.2,3,4,0))
               layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                     ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Find the plot limit for the y scale.                                  #
               #---------------------------------------------------------------------------#
               if (ed22.ci){
                  yuse   = c(sta.q025,sta.q975,ed2.q025,ed2.q975)
                  ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)



                  #----- Plot legend. -----------------------------------------------------#
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
                  #------------------------------------------------------------------------#
               }else{
                  yuse   = c(sta.q025,sta.q975,ed2.median)
                  ylimit = pretty.xylim(u=yuse,fracexp=0.0,is.log=ylog)



                  #----- Plot legend. -----------------------------------------------------#
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
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Loop over all years.                                                 #
               #---------------------------------------------------------------------------#
               for (y in 2:n.census){
                  k       = y - 1
                  left    = (k %% lo.box$ncol) == 1
                  right   = (k %% lo.box$ncol) == 0
                  top     = k <= lo.box$ncol
                  bottom  = k > (lo.box$nrow - 1) * lo.box$ncol
                  mar.now = c(1 + 3 * bottom,1 + 1 * left,1 + 2 * top,1 + 1 * right) + 0.1



                  #------------------------------------------------------------------------#
                  #      95% Confidence Interval.                                          #
                  #------------------------------------------------------------------------#
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
                  #------------------------------------------------------------------------#



                  #----- Set up the title for each plot. ----------------------------------#
                  lesub = paste("Census period: ",census.desc[y],sep="")
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Go on and plot stuff.                                              #
                  #------------------------------------------------------------------------#
                  #----- Plotting window and grid. ----------------------------------------#
                  par(mar=mar.now)
                  plot.new()
                  plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                  if (bottom) axis(side=1,at=x.dbh,labels=dbh.names)
                  if (left  ) axis(side=2)
                  box()
                  title(main=lesub)
                  if (plotgrid) abline(v=x.edge,h=axTicks(2),col=grid.colour,lty="solid")
                  #----- Plot the taxon rate with confidence interval. --------------------#
                  epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col,angle=c(-45,45)
                          ,density=40,lty="solid",lwd=1.0)
                  lines(x=x.dbh,y=sta.median[,y],type="o",col=col.sta[1],pch=16,lwd=2.0)
                  lines(x=x.dbh,y=ed2.median[,y],type="o",col=col.ed2[1],pch=16,lwd=2.0)
               }#end for
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Make the title and axis labels.                                       #
               #---------------------------------------------------------------------------#
               letitre = paste("Size-dependent ",desc.rate,"\n",lieu,sep="")
               ley     = paste(desc.rate," [% ",desc.indiv[i],"/yr]",sep="")
               lex     = "DBH class [cm]"
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot the global title.                                                #
               #---------------------------------------------------------------------------#
               par(las=0)
               mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
               mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
               mtext(text=letitre,side=3,outer=TRUE,cex=cex.main,font=2)
               par(par.orig)
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
            #------------------------------------------------------------------------------#

         }#end for (i in 1:nindiv)
         #---------------------------------------------------------------------------------#
      }#end for
      #====================================================================================#
      #====================================================================================#





      #====================================================================================#
      #====================================================================================#
      #  7. Plot the theme time series of the rates (by DBH class if applicable).          #
      #====================================================================================#
      #====================================================================================#
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
         #---------------------------------------------------------------------------------#
         #    Loop over all possible types of population count.                            #
         #---------------------------------------------------------------------------------#
         for (i in 1:nindiv){
            #==============================================================================#
            #==============================================================================#
            #     PLOT-LEVEL rates.                                                        #
            #------------------------------------------------------------------------------#
            cat("  - Plot level: ",desc.indiv[i],"...","\n")
            ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
            sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


            #----- Create a directory for this type of plot. ------------------------------#
            outtheme   = paste(outplot,theme.now,sep="/")
            if (! file.exists(outtheme)) dir.create(outtheme)
            #------------------------------------------------------------------------------#


            #----- Create path for this individual. ---------------------------------------#
            outindiv   = paste(outtheme,indiv[i],sep="/")
            if (! file.exists(outindiv)) dir.create(outindiv)
            #------------------------------------------------------------------------------#


            #----- Load the modelled rates. -----------------------------------------------#
            sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
            sta.median   = 100. * sta.mod[4,2:n.census]
            sta.q025     = 100. * sta.mod[5,2:n.census]
            sta.q975     = 100. * sta.mod[6,2:n.census]
            ed2.median   = list()
            ed2.q025     = list()
            ed2.q975     = list()
            for (r in sequence(nrate)){
               ed2.mod         = ed2[[ed2.rate[r]]][[indiv[i]]]$global
               ed2.median[[r]] = 100. * ed2.mod[4,2:n.census]
               ed2.q025  [[r]] = 100. * ed2.mod[5,2:n.census]
               ed2.q975  [[r]] = 100. * ed2.mod[6,2:n.census]
            }#end for
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Find the DBH for x scale.                                                 #
            #------------------------------------------------------------------------------#
            x.years  = year4[2:n.census]
            xlimit   = range(x.years)
            #------------------------------------------------------------------------------#


            if (ed22.ci){
               #---------------------------------------------------------------------------#
               #    Make the polygons.                                                     #
               #---------------------------------------------------------------------------#
               plot.poly     = list()
               plot.poly$x   = c(x.years ,rev(x.years) )
               plot.poly$y   = c(sta.q025,rev(sta.q975))
               plot.poly$col = c(col.sta[2])
               for (r in sequence(nrate)){
                  plot.poly$x   = c(plot.poly$x  ,NA,x.years      ,rev(x.years)      )
                  plot.poly$y   = c(plot.poly$y  ,NA,ed2.q025[[r]],rev(ed2.q975[[r]]))
                  plot.poly$col = c(plot.poly$col,col.ed2[r,2])
               }#end for
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find the plot limit for the y scale.                                  #
               #---------------------------------------------------------------------------#
               yuse   = c(sta.q025,sta.q975,unlist(ed2.q025),unlist(ed2.q975))
               ylimit = pretty.xylim(u=yuse,fracexp=0.,is.log=ylog)
               #---------------------------------------------------------------------------#
            }else{
               #---------------------------------------------------------------------------#
               #    Make the polygons.                                                     #
               #---------------------------------------------------------------------------#
               plot.poly     = list()
               plot.poly$x   = c(x.years ,rev(x.years) )
               plot.poly$y   = c(sta.q025,rev(sta.q975))
               plot.poly$col = col.sta[2]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find the plot limit for the y scale.                                  #
               #---------------------------------------------------------------------------#
               yuse   = c(sta.q025,sta.q975,unlist(ed2.median))
               ylimit = pretty.xylim(u=yuse,fracexp=0.,is.log=ylog)
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all formats, and make the plots.                               #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open the file or the plot window. -----------------------------------#
               fichier = paste(outindiv,"/theme-",theme.now,"-",indiv[i],"."
                              ,outform[o],sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Make the title and axis labels.                                       #
               #---------------------------------------------------------------------------#
               letitre = paste(theme.desc," - ",lieu,sep="")
               ley     = paste(theme.desc," [% ",desc.indiv[i],"/yr]",sep="")
               lex     = "Census year"
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Go on and plot stuff.                                                 #
               #---------------------------------------------------------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(5,1))
               #---------------------------------------------------------------------------#


               #----- Plot legend. --------------------------------------------------------#
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
               #---------------------------------------------------------------------------#


               #----- Plotting window and grid. -------------------------------------------#
               par(mar=c(5,4,4,2)+0.1)
               plot(x=x.years,y=sta.median,xlim=xlimit,ylim=ylimit,type="n",main=letitre
                   ,xlab=lex,ylab=ley,log=plog,cex.main=0.7)
               if (plotgrid) grid(col=grid.colour,lty="solid")
               #----- Plot the taxon rate with confidence interval. -----------------------#
               epolygon(x=plot.poly$x,y=plot.poly$y,col=plot.poly$col,angle=c(90,angle)
                       ,density=c(40,dens),lty="solid",lwd=1.0)
               lines(x=x.years,y=sta.median,type="o",col=col.sta[1],pch=16,lwd=2.0)
               for (r in 1:nrate){
                  lines(x=x.years,y=ed2.median[[r]],type="o",col=col.ed2[r,1]
                       ,pch=16,lwd=2.0)
               }#end for
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
            #------------------------------------------------------------------------------#




            #==============================================================================#
            #==============================================================================#
            #     DBH-LEVEL rates.                                                         #
            #------------------------------------------------------------------------------#
            if (sizetoo){
               cat("  - DBH classes: ",desc.indiv[i],"...","\n")
               ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
               sta.rate = paste(this.plot$sta.rate,"size",sep=".")


               #----- Create a directory for this type of plot. ---------------------------#
               outtheme   = paste(outsize,theme.now,sep="/")
               if (! file.exists(outtheme)) dir.create(outtheme)
               #---------------------------------------------------------------------------#


               #----- Create path for this individual. ------------------------------------#
               outindiv   = paste(outtheme,indiv[i],sep="/")
               if (! file.exists(outindiv)) dir.create(outindiv)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    Find the DBH for x scale.                                              #
               #---------------------------------------------------------------------------#
               x.years  = year4[2:n.census]
               xlimit   = range(x.years)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Define a nice configuration for the multiple panels.                 #
               #---------------------------------------------------------------------------#
               lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over all formats, and make the plots.                            #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Open the file or the plot window. --------------------------------#
                  fichier = paste(outindiv,"/theme-",theme.now,"-",indiv[i],".",outform[o]
                                 ,sep="")
                  if(outform[o] == "x11"){
                     X11(width=wide.size$width,height=wide.size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wide.size$width*depth
                        ,height=wide.size$height*depth,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wide.size$width,height=wide.size$height
                               ,pointsize=ptsz,paper=wide.size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE
                        ,width=wide.size$width,height=wide.size$height,pointsize=ptsz
                        ,paper=wide.size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.                     #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  par(oma = c(0.2,3,4,0))
                  layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                        ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
                  par(mar=c(0.1,0.1,0.1,0.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
                  if (ed22.ci){
                     legend ( x       = "bottom"
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
                            , cex     = 0.85
                            , xpd     = TRUE
                            )#end legend
                  }else{
                     legend ( x       = "bottom"
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
                            , xpd     = TRUE
                            )#end legend
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #    Loop over all DBH classes.                                          #
                  #------------------------------------------------------------------------#
                  for (d in 1:n.dbh){
                     left    = (d %% lo.box$ncol) == 1
                     right   = (d %% lo.box$ncol) == 0
                     top     = d <= lo.box$ncol
                     bottom  = d > (lo.box$nrow - 1) * lo.box$ncol




                     #----- Load the modelled rates. --------------------------------------#
                     sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
                     sta.median   = 100. * sta.mod[4,d,2:n.census]
                     sta.q025     = 100. * sta.mod[5,d,2:n.census]
                     sta.q975     = 100. * sta.mod[6,d,2:n.census]
                     ed2.median   = list()
                     ed2.q025     = list()
                     ed2.q975     = list()
                     for (r in 1:nrate){
                        ed2.mod           = ed2[[ed2.rate[r]]][[indiv[i]]]$global
                        ed2.median[[r]]   = 100. * ed2.mod[4,d,2:n.census]
                        ed2.q025  [[r]]   = 100. * ed2.mod[5,d,2:n.census]
                        ed2.q975  [[r]]   = 100. * ed2.mod[6,d,2:n.census]
                     }#end for
                     #---------------------------------------------------------------------#



                     if (ed22.ci){
                        #------------------------------------------------------------------#
                        #     Find the plot limit for the y scale.                         #
                        #------------------------------------------------------------------#
                        yuse   = c(unlist(ed2.q025),unlist(ed2.q975),sta.q025,sta.q975)
                        ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #    Make the polygons.                                            #
                        #------------------------------------------------------------------#
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
                        #------------------------------------------------------------------#
                     }else{
                        #------------------------------------------------------------------#
                        #     Find the plot limit for the y scale.                         #
                        #------------------------------------------------------------------#
                        yuse   = c(ed2.median[d,],sta.q025,sta.q975)
                        ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #    Make the polygons.                                            #
                        #------------------------------------------------------------------#
                        size.poly     = list()
                        size.poly$x   = c(x.years   ,rev(x.years)     )
                        size.poly$y   = c(sta.q025  ,rev(sta.q975)    )
                        size.poly$col = c(col.sta[2])
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#


                     #----- Set up the title and axes labels. -----------------------------#
                     lesub = paste("DBH class:",dbh.names[d],sep="")
                     #---------------------------------------------------------------------#


                     #----- Plot the box plot. --------------------------------------------#
                     par(mar=c(2,2,4,1)+0.1)
                     #----- Plotting window and grid. -------------------------------------#
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                     axis(side=1)
                     axis(side=2)
                     box()
                     title(main=lesub,xlab="",ylab="")
                     if (plotgrid) grid(col=grid.colour,lty="solid")
                     #----- Plot the taxon rate with confidence interval. -----------------#
                     epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col
                             ,angle=c(90,angle),density=c(40,dens),lty="solid",lwd=1.0)
                     lines(x=x.years,y=sta.median,type="o",pch=16,lwd=2.0
                          ,col=col.sta[1])
                     for (r in sequence(nrate)){
                        lines(x=x.years,y=ed2.median[[r]],type="o",pch=16,lwd=2.0
                             ,col=col.ed2[r,1])
                     }#end for
                     #---------------------------------------------------------------------#
                  }#end for (d in 1:n.dbh)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Make the title and axis labels.                                    #
                  #------------------------------------------------------------------------#
                  letitre = paste(theme.desc,": ",lieu,sep="")
                  ley     = paste(theme.desc," [% ",desc.indiv[i],"]",sep="")
                  lex     = "Census"
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the plotting window.                                         #
                  #------------------------------------------------------------------------#
                  par(las=0)
                  mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  mtext(text=letitre,side=3,outer=TRUE,cex=cex.main,font=2)
                  par(par.orig)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #====================================================================================#
      #====================================================================================#





      #====================================================================================#
      #====================================================================================#
      #  8. Plot the time series of the rates (by DBH class if applicable).                #
      #====================================================================================#
      #====================================================================================#
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
         #---------------------------------------------------------------------------------#
         #    Loop over all possible types of population count.                            #
         #---------------------------------------------------------------------------------#
         for (i in 1:nindiv){
            #==============================================================================#
            #==============================================================================#
            #     PLOT-LEVEL rates.                                                        #
            #------------------------------------------------------------------------------#
            cat("  - Plot level: ",desc.indiv[i],"...","\n")
            ed2.rate       = paste(this.plot$ed2.rate,"plot",sep=".")
            sta.rate       = paste(this.plot$sta.rate,"plot",sep=".")


            #----- Create a directory for this type of plot. ------------------------------#
            outrate   = paste(outplot,ed2.rate,sep="/")
            if (! file.exists(outrate)) dir.create(outrate)
            #------------------------------------------------------------------------------#


            #----- Create path for this individual. ---------------------------------------#
            outindiv   = paste(outrate,indiv[i],sep="/")
            if (! file.exists(outindiv)) dir.create(outindiv)
            #------------------------------------------------------------------------------#


            #----- Load the modelled rates. -----------------------------------------------#
            sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
            sta.median   = 100. * sta.mod[4,2:n.census]
            sta.q025     = 100. * sta.mod[5,2:n.census]
            sta.q975     = 100. * sta.mod[6,2:n.census]
            ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
            ed2.median   = 100. * ed2.mod[4,2:n.census]
            ed2.q025     = 100. * ed2.mod[5,2:n.census]
            ed2.q975     = 100. * ed2.mod[6,2:n.census]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Find the DBH for x scale.                                                 #
            #------------------------------------------------------------------------------#
            x.years  = year4[2:n.census]
            xlimit   = range(x.years)
            #------------------------------------------------------------------------------#


            if (ed22.ci){
               #---------------------------------------------------------------------------#
               #    Make the polygons.                                                     #
               #---------------------------------------------------------------------------#
               plot.poly     = list()
               plot.poly$x   = c(x.years ,rev(x.years) ,NA,x.years ,rev(x.years) )
               plot.poly$y   = c(sta.q025,rev(sta.q975),NA,ed2.q025,rev(ed2.q975))
               plot.poly$col = c(col.sta[2],col.ed2[2])
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find the plot limit for the y scale.                                  #
               #---------------------------------------------------------------------------#
               yuse   = c(sta.q025,sta.q975,ed2.q025,ed2.q975)
               ylimit = pretty.xylim(u=yuse,fracexp=scalleg,is.log=ylog)
               #---------------------------------------------------------------------------#
            }else{
               #---------------------------------------------------------------------------#
               #    Make the polygons.                                                     #
               #---------------------------------------------------------------------------#
               plot.poly     = list()
               plot.poly$x   = c(x.years ,rev(x.years) )
               plot.poly$y   = c(sta.q025,rev(sta.q975))
               plot.poly$col = col.sta[2]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Find the plot limit for the y scale.                                  #
               #---------------------------------------------------------------------------#
               yuse   = c(sta.q025,sta.q975,ed2.median)
               ylimit = pretty.xylim(u=yuse,fracexp=scalleg,is.log=ylog)
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over all formats, and make the plots.                               #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open the file or the plot window. -----------------------------------#
               fichier = paste(outindiv,"/tseries-",ed2.rate,"-",indiv[i],"."
                              ,outform[o],sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#
               
               par(par.user)
               layout(mat=rbind(2,1),heights=c(5,1))



               #----- Plot legend. --------------------------------------------------------#
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
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Make the title and axis labels.                                       #
               #---------------------------------------------------------------------------#
               letitre = paste(desc.rate," - ",lieu,sep="")
               ley     = paste(desc.rate," [% ",desc.indiv[i],"/yr]",sep="")
               lex     = "Census year"
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Go on and plot stuff.                                                 #
               #---------------------------------------------------------------------------#
               par(mar=c(5,4,4,2)+0.1)
               #----- Plotting window and grid. -------------------------------------------#
               plot(x=x.years,y=sta.median,xlim=xlimit,ylim=ylimit,type="n",main=letitre
                   ,xlab=lex,ylab=ley,log=plog,cex.main=0.7)
               if (plotgrid) grid(col=grid.colour,lty="solid")
               #----- Plot the taxon rate with confidence interval. -----------------------#
               epolygon(x=plot.poly$x,y=plot.poly$y,col=plot.poly$col,angle=c(-45,45)
                       ,density=40,lty="solid",lwd=1.0)
               lines(x=x.years,y=sta.median,type="o",col=col.sta[1],pch=16,lwd=2.0)
               lines(x=x.years,y=ed2.median,type="o",col=col.ed2[1],pch=16,lwd=2.0)


               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#




            #==============================================================================#
            #==============================================================================#
            #     DBH-LEVEL rates.                                                         #
            #------------------------------------------------------------------------------#
            if (sizetoo){
               cat("  - DBH classes: ",desc.indiv[i],"...","\n")
               ed2.rate = paste(this.plot$ed2.rate,"size",sep=".")
               sta.rate = paste(this.plot$sta.rate,"size",sep=".")


               #----- Create a directory for this type of plot. ---------------------------#
               outrate   = paste(outsize,ed2.rate,sep="/")
               if (! file.exists(outrate)) dir.create(outrate)
               #---------------------------------------------------------------------------#


               #----- Create path for this individual. ------------------------------------#
               outindiv   = paste(outrate,indiv[i],sep="/")
               if (! file.exists(outindiv)) dir.create(outindiv)
               #---------------------------------------------------------------------------#


               #----- Load the modelled rates. --------------------------------------------#
               sta.mod      = sta[[sta.rate]][[indiv[i]]]$global
               sta.median   = 100. * sta.mod[4,,2:n.census]
               sta.q025     = 100. * sta.mod[5,,2:n.census]
               sta.q975     = 100. * sta.mod[6,,2:n.census]
               ed2.mod      = ed2[[ed2.rate]][[indiv[i]]]$global
               ed2.median   = 100. * ed2.mod[4,,2:n.census]
               ed2.q025     = 100. * ed2.mod[5,,2:n.census]
               ed2.q975     = 100. * ed2.mod[6,,2:n.census]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    Find the DBH for x scale.                                              #
               #---------------------------------------------------------------------------#
               x.years  = year4[2:n.census]
               xlimit   = range(x.years)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Define a nice configuration for the multiple panels.                 #
               #---------------------------------------------------------------------------#
               lo.box = pretty.box(n=n.dbh,horizontal=TRUE)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop over all formats, and make the plots.                            #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Open the file or the plot window. --------------------------------#
                  fichier = paste(outindiv,"/tseries-",ed2.rate,"-",indiv[i],".",outform[o]
                                 ,sep="")
                  if(outform[o] == "x11"){
                     X11(width=wide.size$width,height=wide.size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=wide.size$width*depth
                        ,height=wide.size$height*depth,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=wide.size$width,height=wide.size$height
                               ,pointsize=ptsz,paper=wide.size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE
                        ,width=wide.size$width,height=wide.size$height,pointsize=ptsz
                        ,paper=wide.size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the window into several smaller windows.                     #
                  #------------------------------------------------------------------------#
                  par(par.user)
                  par.orig = par(no.readonly = TRUE)
                  par(oma = c(0.2,3,4,0))
                  layout(mat    = rbind(1+lo.box$mat,rep(1,times=lo.box$ncol))
                        ,height = c(rep(5/lo.box$nrow,times=lo.box$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#



                  #----- Plot legend. -----------------------------------------------------#
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
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #    Loop over all DBH classes.                                          #
                  #------------------------------------------------------------------------#
                  for (d in 1:n.dbh){
                     left    = (d %% lo.box$ncol) == 1
                     right   = (d %% lo.box$ncol) == 0
                     top     = d <= lo.box$ncol
                     bottom  = d > (lo.box$nrow - 1) * lo.box$ncol

                     if (ed22.ci){
                        #------------------------------------------------------------------#
                        #     Find the plot limit for the y scale.                         #
                        #------------------------------------------------------------------#
                        yuse   = c(ed2.q025[d,],ed2.q975[d,],sta.q025[d,],sta.q975[d,])
                        ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #    Make the polygons.                                            #
                        #------------------------------------------------------------------#
                        size.poly  = list()
                        size.poly$x   = c(x.years     ,rev(x.years)     ,NA
                                         ,x.years     ,rev(x.years)         )
                        size.poly$y   = c(sta.q025[d,],rev(sta.q975[d,]),NA
                                         ,ed2.q025[d,],rev(ed2.q975[d,])    )
                        size.poly$col = c(col.sta[2]  ,col.ed2[2]           )
                        #------------------------------------------------------------------#
                     }else{
                        #------------------------------------------------------------------#
                        #     Find the plot limit for the y scale.                         #
                        #------------------------------------------------------------------#
                        yuse   = c(ed2.median[d,],sta.q025[d,],sta.q975[d,])
                        ylimit = pretty.xylim(yuse,fracexp=0.0,is.log=ylog)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #    Make the polygons.                                            #
                        #------------------------------------------------------------------#
                        size.poly     = list()
                        size.poly$x   = c(x.years     ,rev(x.years)         )
                        size.poly$y   = c(sta.q025[d,],rev(sta.q975[d,])    )
                        size.poly$col = c(col.sta[2]  )
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#


                     #----- Set up the title and axes labels. -----------------------------#
                     lesub = paste("DBH class:",dbh.names[d],sep="")
                     #---------------------------------------------------------------------#


                     #----- Plot the box plot. --------------------------------------------#
                     par(mar=c(2,2,4,1)+0.1)
                     #----- Plotting window and grid. -------------------------------------#
                     plot.new()
                     plot.window(xlim=xlimit,ylim=ylimit,log=plog)
                     axis(side=1)
                     axis(side=2)
                     box()
                     title(main=lesub,xlab="",ylab="")
                     if (plotgrid) grid(col=grid.colour,lty="solid")
                     #----- Plot the taxon rate with confidence interval. -----------------#
                     epolygon(x=size.poly$x,y=size.poly$y,col=size.poly$col,angle=c(-45,45)
                             ,density=40,lty="solid",lwd=1.0)
                     lines(x=x.years,y=sta.median[d,],type="o",pch=16,lwd=2.0
                          ,col=col.sta[1])
                     lines(x=x.years,y=ed2.median[d,],type="o",pch=16,lwd=2.0
                          ,col=col.ed2[1])
                     #---------------------------------------------------------------------#
                  }#end for (d in 1:n.dbh)
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Make the title and axis labels.                                    #
                  #------------------------------------------------------------------------#
                  letitre = paste(desc.rate,": ",lieu,sep="")
                  ley     = paste(desc.rate," [% ",desc.indiv[i],"]",sep="")
                  lex     = "Census"
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Split the plotting window.                                         #
                  #------------------------------------------------------------------------#
                  par(las=0)
                  mtext(text=lex    ,side=1,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
                  mtext(text=ley    ,side=2,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
                  mtext(text=letitre,side=3,outer=TRUE,cex=cex.main,font=2)
                  par(par.orig)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #====================================================================================#
      #====================================================================================#
   }#end if (census.name %in% ls())
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #=======================================================================================#
   #     Plot the monthly mean variables as functions of other 2 environment variables.    #
   #---------------------------------------------------------------------------------------#
   cat(" + Plotting parameter space...",sep="","\n")


   #----- Sizes. --------------------------------------------------------------------------#
   nzvar = length(xyzvar$zvar)
   nxvar = length(xyzvar$xvar)
   nyvar = length(xyzvar$yvar)
   #---------------------------------------------------------------------------------------#


   #----- Create the directories. ---------------------------------------------------------#
   outxyzp = paste(outpref,"xyzplot",sep="/")
   if (! file.exists(outxyzp )) dir.create(outxyzp )
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    Loop over all explanatory variables that go to the Y axis.                         #
   #---------------------------------------------------------------------------------------#
   for (y in 1:nyvar){
      #----- Load Y settings. -------------------------------------------------------------#
      this.y  = xyzvar$yvar[[y]]
      yvname  = this.y$vname
      ydesc   = this.y$desc
      yunit   = this.y$unit
      yadd    = this.y$add
      ymult   = this.y$mult
      yleg    = this.y$leg
      ylog    = this.y$log
      sizetoo = this.y$sizetoo
      cat("     * Y: ",ydesc,"...","\n")
      #------------------------------------------------------------------------------------#



      #----- Create the directories. ------------------------------------------------------#
      outyvar = paste(outxyzp ,yvname   ,sep="/")
      if (! file.exists(outyvar )) dir.create(outyvar )
      #------------------------------------------------------------------------------------#



      #----- Load the y variable.  Expand the edges of the y axis to fit a legend. --------#
      if (sizetoo){
         yplot.name  = paste(yvname,"plot","mean",sep=".")
         ysize.name  = paste(yvname,"size","mean",sep=".")
      }else{
         yplot.name  = paste(yvname,"mean",sep=".")
         ysize.name  = paste(yvname,"mean",sep=".")
      }#end if
      yvar.plot   = ymult * ( pspace[[yplot.name]] + yadd )
      yvar.size   = ymult * ( pspace[[ysize.name]] + yadd )
      ley         = paste(ydesc," [",yunit,"]",sep="")
      ylimit.plot = pretty.xylim(u=yvar.plot,fracexp=scalleg,is.log=ylog)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Loop over all explanatory variables that go to the X axis.                      #
      #------------------------------------------------------------------------------------#
      for (x in 1:nxvar){
         #----- Load X settings. ----------------------------------------------------------#
         this.x = xyzvar$xvar[[x]]
         xvname = this.x$vname
         xdesc  = this.x$desc
         xunit  = this.x$unit
         xadd   = this.x$add
         xmult  = this.x$mult
         xleg   = this.x$leg
         xlog   = this.x$log
         #---------------------------------------------------------------------------------#



         #----- Load the x variable. ------------------------------------------------------#
         xplot.name  = paste(xvname,"plot","mean",sep=".")
         xsize.name  = paste(xvname,"size","mean",sep=".")
         xvar.plot   = xmult * ( pspace[[xplot.name]] + xadd )
         xvar.size   = xmult * ( pspace[[xsize.name]] + xadd )
         lex         = paste(xdesc," [",xunit,"]",sep="")

         sel         = is.finite(xvar.plot) & ( xvar.plot > 0 | (! xlog))
         xlimit.plot = range(xvar.plot[sel],na.rm=TRUE)
         xlimit.size = matrix(nrow=n.dbh,ncol=2)
         for (d in 1:n.dbh){
            sel             = is.finite(xvar.size[d,]) & ( xvar.size[d,] > 0 | (! xlog) )
            xlimit.size[d,] = range(xvar.size[d,sel],na.rm=TRUE)
         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #   Make the log scale, and find the position for the legend.                     #
         #---------------------------------------------------------------------------------#
         plog = ""
         if (xlog) plog=paste(plog,"x",sep="")
         if (ylog) plog=paste(plog,"y",sep="")
         leg.pos   = paste(yleg,xleg,sep="")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Loop over all colour variables.                                              #
         #---------------------------------------------------------------------------------#
         for (z in 1:nzvar){
            #----- Load Z settings. -------------------------------------------------------#
            this.z   = xyzvar$zvar[[z]]
            zvname   = this.z$vname
            zdesc    = this.z$desc
            zkey     = this.z$key
            zunit    = this.z$unit
            zadd     = this.z$add
            zmult    = this.z$mult
            zcscheme = get(this.z$col.scheme)
            zlog     = this.z$log
            cat("      ~ X:",xdesc,"   Z: ",zdesc,"...","\n")
            #------------------------------------------------------------------------------#



            #----- Annotation for the colour map ("Z" axis). ------------------------------#
            zvar  = zmult * ( pspace[[zvname]] + zadd )
            lez   = paste(zkey,"\n [",zunit,"]",sep="")
            #------------------------------------------------------------------------------#



            #----- Find the range for the scale. ------------------------------------------#
            zlimit  = pretty.xylim(u=zvar,fracexp=0.0,is.log=zlog)
            #------------------------------------------------------------------------------#


            #----- Title. -----------------------------------------------------------------#
            letitre.plot = paste(lieu,paste(zdesc,"Plot level",sep=" - "),sep="\n")
            letitre.size = paste(lieu,paste(zdesc,"Size level",sep=" - "),sep="\n")
            lesub.size   = paste("DBH Class: ",dbh.names,sep="")
            #------------------------------------------------------------------------------#



            #----- Attribute symbols according to the year. -------------------------------#
            this.pch  = pspace$pch
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the plot-level.                                                     #
            #------------------------------------------------------------------------------#
            cat("         > Plot-level...","\n")
            for (o in 1:nout){
               #----- Open the file. ------------------------------------------------------#
               fichier = paste(outyvar,"/plot_x_",xvname,"_y_",yvname
                                      ,"_z_",zvname,".",outform[o],sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#



               #----- Plot the parameter space. -------------------------------------------#
               par(par.user)
               xyz.plot( x              = xvar.plot
                       , y              = yvar.plot
                       , z              = zvar
                       , fixed.xlim     = TRUE
                       , fixed.ylim     = TRUE
                       , xlim           = xlimit.plot
                       , ylim           = ylimit.plot
                       , zlim           = zlimit
                       , colour.palette = zcscheme
                       , cex            = 1.6
                       , pch            = this.pch
                       , xy.log         = plog
                       , xyz.main       = list(text=letitre.plot,cex=cex.main)
                       , xyz.xlab       = list(text=lex,adj=mtext.xadj,padj=mtext.xoff)
                       , xyz.ylab       = list(text=ley,adj=mtext.yadj,padj=mtext.yoff)
                       , key.title      = list(main=lez,cex.main=0.8)
                       , key.log        = zlog
                       , xyz.more       = list(grid=list(col=grid.colour,lty="solid"))
                       , xyz.legend     = list( x      = "bottom"
                                              , inset  = 0.01
                                              , legend = pspace$leg.label
                                              , col    = foreground
                                              , bg     = background
                                              , pch    = pspace$leg.pch
                                              , title  = "Census"
                                              , ncol   = 3
                                              , cex    = 1.0
                                              )#end legend
                       )#end xyz.plot
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Close the device.                                                     #
               #---------------------------------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the parameter space for the dbh classes if it is supposed to do so. #
            #------------------------------------------------------------------------------#
            if (sizetoo){
               xuse.size = split(x=xvar.size,f=col(xvar.size))
               yuse.size = split(x=yvar.size,f=col(yvar.size))
               pch.size  = replicate(n=length(xuse.size),expr=c(this.pch),simplify=FALSE)


               #---------------------------------------------------------------------------#
               #     Plot the bar plot.                                                    #
               #---------------------------------------------------------------------------#
               for (o in 1:nout){
                  #----- Open the file. ---------------------------------------------------#
                  fichier = paste(outyvar,"/size_x_",xvname
                                 ,"_y_",yvname,"_z_",zvname,".",outform[o],sep="")
                  if(outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth
                        ,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=size$paper)
                  }else if(outform[o] == "pdf"){
                     pdf(file=fichier,onefile=FALSE
                        ,width=size$width,height=size$height,pointsize=ptsz
                        ,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Plot the parameter space. ----------------------------------------#
                  par(par.user)
                  xyz.plot( x              = xuse.size
                          , y              = yuse.size
                          , z              = zvar
                          , fixed.xlim     = FALSE
                          , fixed.ylim     = FALSE
                          , zlim           = zlimit
                          , colour.palette = muitas
                          , cex            = 1.6
                          , pch            = this.pch
                          , xy.log         = plog
                          , xyz.main       = list(text=letitre.size,cex.main=cex.main)
                          , xyz.xlab       = list(text=lex,adj=mtext.xadj,padj=mtext.xoff)
                          , xyz.ylab       = list(text=ley,adj=mtext.yadj,padj=mtext.yoff)
                          , key.title      = list(main=lez,cex.main=0.8)
                          , key.log        = zlog
                          , xyz.more       = list(grid=list(col=grid.colour,lty="solid"))
                          , xyz.legend     = list( x      = "bottom"
                                                 , inset  = 0.01
                                                 , legend = pspace$leg.label
                                                 , col    = foreground
                                                 , bg     = background
                                                 , pch    = pspace$leg.pch
                                                 , title  = "Census"
                                                 , ncol   = 3
                                                 , cex    = 1.0
                                                 )#end legend
                          )#end xyz.plot
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Close the device.                                                  #
                  #------------------------------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for (o in 1:nout)
               #---------------------------------------------------------------------------#
            }#end if (sizetoo)
            #------------------------------------------------------------------------------#
         }#end for (x in 1:nxvar)
         #---------------------------------------------------------------------------------#
      }#end for (y in 1:nyvar)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#
