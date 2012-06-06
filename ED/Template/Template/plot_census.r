rm(list=ls())

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
here           = "thispath"    # Current directory.
there          = "thatpath"    # Directory where analyses/history are 
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc"      # Source  directory.
outroot        = "thisoutroot"
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
various.cycles = myvarcycle
myplaces       = c("thispoly")
sasmonth.short = c(2,5,8,11)
sasmonth.long  = 5
nyears.long    = 25
outform        = "thisoutform"          # Formats for output file.  Supported formats are:
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
                                #   - "pdf" - for PDF printing

byeold         = TRUE           # Remove old files of the given format?

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.20
cex.main       = 0.8            # Scale coefficient for the title
slz.min        = -5.0           # Find the deepest depth that trees access water.
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



#------------------------------------------------------------------------------------------#
#     Plot-level comparisons.                                                              #
#------------------------------------------------------------------------------------------#
plotvar       = list()
plotvar[[ 1]] = list( vnam.ed    = "recr"
                    , vnam.obs   = "recr"
                    , desc       = "Recruitment rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("chartreuse4","olivedrab2")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 2]] = list( vnam.ed    = "mort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Total mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("purple4","mediumpurple1")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 3]] = list( vnam.ed    = "ddmort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Density dependent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("orangered","orange")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 4]] = list( vnam.ed    = "dimort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Density independent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("sienna4","darkgoldenrod2")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 5]] = list( vnam.ed    = "growth.plot"
                    , vnam.obs   = "growth.plot"
                    , desc       = "Growth rate"
                    , unit       = "[%DBH/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("royalblue4","lightskyblue")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Plot-level comparisons.                                                              #
#------------------------------------------------------------------------------------------#
sizevar       = list()
sizevar[[ 1]] = list( vnam.ed    = "mort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Total mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("purple4","mediumpurple1")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 2]] = list( vnam.ed    = "ddmort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Density dependent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("orangered","orange")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 3]] = list( vnam.ed    = "dimort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Density independent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("sienna4","darkgoldenrod2")
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 4]] = list( vnam.ed    = "growth.size"
                    , vnam.obs   = "growth.size"
                    , desc       = "Growth rate"
                    , unit       = "[%DBH/yr]"
                    , col.obser  = c("gray21","gray42")
                    , col.model  = c("royalblue4","lightskyblue")
                    , leg.corner = "topleft"
                    , plog       = TRUE
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
xyzvar$yvar      = list( list( vname   = "recr"
                             , desc    = "Recruitment rate"
                             , key     = "Recruitment"
                             , unit    = "%pop/yr"
                             , add     = 0.
                             , mult    = 100.
                             , leg     = "top"
                             , log     = TRUE
                             , sizetoo = FALSE
                             )#end list
                       , list( vname   = "mort"
                             , desc    = "Mortality rate"
                             , key     = "Mortality"
                             , unit    = "%pop/yr"
                             , add     = 0.
                             , mult    = 100.
                             , leg     = "top"
                             , log     = TRUE
                             , sizetoo = TRUE
                             )#end list
                       , list( vname   = "ddmort"
                             , desc    = "Density-dependent mortality rate"
                             , key     = "DD Mort."
                             , unit    = "%pop/yr"
                             , add     = 0.
                             , mult    = 100.
                             , leg     = "top"
                             , log     = TRUE
                             , sizetoo = TRUE
                             )#end list
                       , list( vname   = "dimort"
                             , desc    = "Density-independent mortality rate"
                             , key     = "DI Mort."
                             , unit    = "%pop/yr"
                             , add     = 0.
                             , mult    = 100.
                             , leg     = "top"
                             , log     = TRUE
                             , sizetoo = TRUE
                             )#end list
                       , list( vname   = "growth"
                             , desc    = "Growth rate"
                             , key     = "Growth"
                             , unit    = "%pop/yr"
                             , add     = 0.
                             , mult    = 100.
                             , leg     = "top"
                             , log     = TRUE
                             , sizetoo = TRUE
                             )#end list
                       )#end list
xyzvar$zvar      = list( list( vname = "rshort"
                             , desc  = "Mean shortwave radiation"
                             , unit  = "W/m2"
                             , add   = 0.
                             , mult  = 1.
                             , log   = FALSE
                             )#end list
                       , list( vname = "fs.open"
                             , desc  = "Minimum water stress scale"
                             , unit  = "--"
                             , add   = 0.
                             , mult  = 1.
                             , log   = FALSE
                             )#end list
                       , list( vname = "paw"
                             , desc  = "Minimum available water"
                             , unit  = "%"
                             , add   = 0.
                             , mult  = 100.
                             , log   = FALSE
                             )#end list
                       , list( vname = "smpot"
                             , desc  = "Maximum matric potential"
                             , unit  = "MPa"
                             , add   = 0.
                             , mult  = 1.
                             , log   = FALSE
                             )#end list
                       )#end if


#----- Load some packages. ----------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(Hmisc)
library(sn)
#------------------------------------------------------------------------------------------#


#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#


#----- Set how many variables we will compare. --------------------------------------------#
nplotvar = length(plotvar)
nsizevar = length(sizevar)
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#


#----- Load some files with functions. ----------------------------------------------------#
source(paste(srcdir,"atlas.r"           ,sep="/"))
source(paste(srcdir,"charutils.r"       ,sep="/"))
source(paste(srcdir,"census.r"          ,sep="/"))
source(paste(srcdir,"colourmap.r"       ,sep="/"))
source(paste(srcdir,"cloudy.r"          ,sep="/"))
source(paste(srcdir,"epolygon.r"        ,sep="/"))
source(paste(srcdir,"error.bar.r"       ,sep="/"))
source(paste(srcdir,"globdims.r"        ,sep="/"))
source(paste(srcdir,"locations.r"       ,sep="/"))
source(paste(srcdir,"muitas.r"          ,sep="/"))
source(paste(srcdir,"plotsize.r"        ,sep="/"))
source(paste(srcdir,"pretty.log.r"      ,sep="/"))
source(paste(srcdir,"pretty.time.r"     ,sep="/"))
source(paste(srcdir,"qapply.r"          ,sep="/"))
source(paste(srcdir,"rconstants.r"      ,sep="/"))
source(paste(srcdir,"skewnorm.stats.r"  ,sep="/"))
source(paste(srcdir,"soilutils.r"       ,sep="/"))
source(paste(srcdir,"sombreado.r"       ,sep="/"))
source(paste(srcdir,"southammap.r"      ,sep="/"))
source(paste(srcdir,"thermlib.r"        ,sep="/"))
source(paste(srcdir,"timeutils.r"       ,sep="/"))
#----- These should be called after the others. -------------------------------------------#
source(paste(srcdir,"pft.coms.r"        ,sep="/"))
#------------------------------------------------------------------------------------------#



#----- Load census data. ------------------------------------------------------------------#
census.file = paste(srcdir,"LBA_MIP.census_summ.RData",sep="/")
load(file=census.file)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
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


      #----- Initialise the model structure. ----------------------------------------------#
      model         = list()
      #------------------------------------------------------------------------------------#


      #----- Initialise the parameter space structure. ------------------------------------#
      pspace        = list()
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Load the census data, from the monthly means.                                  #
      #------------------------------------------------------------------------------------#
      if (! "recr" %in% names(model)){
         print(paste("   - Starting new census assessment...",sep=""))
      }else{
         print(paste("   - Resuming census assessment...",sep=""))
      }#end if
      census.obs = get(census.name)
      n.census   = length(census.obs$when)
      dyear      = median(diff(numyears(census.obs$when)))
      n.dbh      = length(census.obs$dbh.breaks)-1
      x.dbh      = c(10,census.obs$dbh.breaks[seq(from=2,to=n.dbh,by=1)])
      dbh.names  = dimnames(census.obs$mort.size)[[2]]


      #------------------------------------------------------------------------------------#
      #      Loop over all months to grab all the census data.                             #
      #------------------------------------------------------------------------------------#
      census.idx   = NULL
      for (y in 2:n.census){
         #----- Find the first and last time to be averaged for this census. --------------#
         ts.montha  = ( nummonths(census.obs$when[y-1]) %% 12 )
         ts.yeara   = numyears (census.obs$when[y-1])
         ts.monthz  = ( ( (nummonths(census.obs$when[y]) - 1) %% 12 )
                      + 12 * as.integer(nummonths(census.obs$when[y]) == 1) )
         ts.yearz   = numyears (census.obs$when[y]) - as.integer(ts.monthz == 12)
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
      act.census.yeara    = numyears(census.obs$when[2])
      act.census.yearz    = numyears(census.obs$when[n.census])
      act.census.year     = seq(from=act.census.yeara,to=act.census.yearz,by=dyear)
      #------------------------------------------------------------------------------------#

      ncyc                = act.census.yearz - act.census.yeara + dyear
      if (various.cycles){
         nfullcyc         = 1 + floor( (act.census.yeara - dyear - yeara ) / ncyc )
      }else{
         nfullcyc         = 1
      }#end if
      i1stfull            = act.census.yeara - (nfullcyc - 1) * ncyc

      n.months = length(census.idx)
      #------------------------------------------------------------------------------------#


      #----- Initialise all the structures for which we will compare. ---------------------#
      ts.recr          = array( NA, dim = c(      n.months,nfullcyc))
      ts.agb.plot      = array( NA, dim = c(      n.months,nfullcyc))
      ts.ba.plot       = array( NA, dim = c(      n.months,nfullcyc))
      ts.lai.plot      = array( NA, dim = c(      n.months,nfullcyc))
      ts.mort.plot     = array( NA, dim = c(      n.months,nfullcyc))
      ts.ddmort.plot   = array( NA, dim = c(      n.months,nfullcyc))
      ts.dimort.plot   = array( NA, dim = c(      n.months,nfullcyc))
      ts.growth.plot   = array( NA, dim = c(      n.months,nfullcyc))
      ts.agb.size      = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      ts.ba.size       = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      ts.lai.size      = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      ts.mort.size     = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      ts.ddmort.size   = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      ts.dimort.size   = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      ts.growth.size   = array( NA, dim = c(n.dbh,n.months,nfullcyc))
      #----- Environmental variables. -----------------------------------------------------#
      ts.rshort        = array( NA, dim = c(      n.months,nfullcyc))
      ts.fs.open       = array( NA, dim = c(      n.months,nfullcyc))
      ts.paw           = array( NA, dim = c(      n.months,nfullcyc))
      ts.smpot         = array( NA, dim = c(      n.months,nfullcyc))
      ts.census.year   = array( NA, dim = c(      n.months,nfullcyc))
      ts.census.idx    = array( NA, dim = c(      n.months,nfullcyc))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop over all times, and retrieve the data.                                    #
      #------------------------------------------------------------------------------------#
      for (u in 1:nfullcyc){
         now.month = nummonths(census.obs$when[1])
         now.year  = i1stfull - dyear + (u-1) * ncyc
         year.use  = act.census.yeara - dyear
         for (m in 1:n.months){
            now.month   = (now.month %% 12) + 1
            now.year    = now.year + as.integer(now.month == 1)
            year.use    = year.use + as.integer(now.month == 1)
            print (paste("ED-2.1: ",paste(mon2mmm(now.month,cap1=T),now.year,sep="-")
                        ,"; Census: ",paste(year.use)))

            #----- Build the file name. ---------------------------------------------------#
            cmonth    = sprintf("%2.2i",now.month)
            cyear     = sprintf("%2.2i",now.year )
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
               print (paste("     * Reading ",basename(myfile),"...",sep=""))
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
               soil.water      = rev(cumsum(rev(mymont$MMEAN.SOIL.WATER * wdns * dslz)))
               soil.moist.avg  = soil.water / (wdns * soil.depth)
               ts.paw    [m,u] = ( ( soil.water[ka] - soil.dry [ka] )
                                 / ( soil.poro [ka] - soil.dry [ka] ) )
               ts.rshort [m,u] = mymont$MMEAN.RSHORT
               ts.fs.open[m,u] = mymont$MMEAN.FS.OPEN
               ts.smpot  [m,u] = ( - smoist2mpot(smoist=soil.moist.avg[ka],mysoil=soil)
                                 * 0.001 * grav )
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
                  dbhcut       = cut(dbhconow,breaks=census.obs$dbh.breaks)
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
                  mortconow       = rowSums(mymont$MMEAN.MORT.RATE)
                  ddmortconow     = mymont$MMEAN.MORT.RATE[,2]
                  dimortconow     = mortconow - ddmortconow
                  recruitconow    = mymont$RECRUIT.DBH
                  censtatusconow  = mymont$CENSUS.STATUS
                  growthconow     = mymont$DLNDBH.DT
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
                  growthconow     = NA
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
               #     Find the growth, mortality, and recruitment rates for each PFT, and   #
               # the global rates.  We only use the cohorts that were flagged as 1 or 2    #
               # (which means that their DBH is greater than 10 cm).                       #
               #---------------------------------------------------------------------------#
               agbpft    = rep(0.,times=npft)
               bapft     = rep(0.,times=npft)
               laipft    = rep(0.,times=npft)
               recrpft   = rep(NA,times=npft)
               mortpft   = rep(NA,times=npft)
               ddmortpft = rep(NA,times=npft)
               dimortpft = rep(NA,times=npft)
               growthpft = rep(NA,times=npft)
               nplantpft = rep(0 ,times=npft)

               for (p in 1:npft){
                  if (all(is.na(pftconow))){
                     sel.dbh      = rep(FALSE,times=length(pftconow))
                     sel.cs2      = rep(FALSE,times=length(pftconow))
                  }else{
                     sel.dbh      = pftconow == p & censtatusconow >  0
                     sel.cs2      = pftconow == p & censtatusconow == 2
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #    Find the AGB, LAI, and BA.                                          #
                  #------------------------------------------------------------------------#
                  if (any(sel.dbh)){
                     agbpft        [p] = sum(w.nplant[sel.dbh] * agbconow[sel.dbh])
                     bapft         [p] = sum(w.nplant[sel.dbh] * baconow [sel.dbh])
                     laipft        [p] = sum(w.lai   [sel.dbh] )
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Find the PFT-level mortality rates.                                #
                  #------------------------------------------------------------------------#
                  if (any(sel.cs2)){

                     #---------------------------------------------------------------------#
                     #    Find the weight of each PFT.                                     #
                     #---------------------------------------------------------------------#
                     nplantpft     [p] = sum( w.nplant[sel.cs2] )
                     #---------------------------------------------------------------------#

                     #---- This is the number of survivors. -------------------------------#
                     survivor          = sum(w.nplant[sel.cs2])
                     previous          = sum(w.nplant[sel.cs2] * exp(mortconow[sel.cs2]))
                     mortpft       [p] = log( previous / survivor )

                     survivor          = sum(w.nplant[sel.cs2])
                     previous          = sum(w.nplant[sel.cs2] * exp(dimortconow[sel.cs2]))
                     dimortpft     [p] = log( previous / survivor )

                     survivor          = sum(w.nplant[sel.cs2])
                     previous          = sum(w.nplant[sel.cs2] * exp(ddmortconow[sel.cs2]))
                     ddmortpft     [p] = log( previous / survivor )
                  }#end if
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Find the PFT-level recruitment rates.                              #
                  #------------------------------------------------------------------------#
                  if (any(sel.dbh) & any(sel.cs2)){
                     #---- This is the number of survivors. -------------------------------#
                     population        = sum(w.nplant[sel.dbh])
                     established       = sum(w.nplant[sel.cs2])
                     recrpft       [p] = log( population / established) / 12.0
                  }#end if
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Growth rates are found only for established cohorts.               #
                  #------------------------------------------------------------------------#
                  if (any(sel.cs2)){
                     growthpft     [p] = sum( w.nplant[sel.cs2] * growthconow [sel.cs2] )
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Find the mean rates amongst all PFTs.  Because Condit et al. (2006)  #
               # assumed log-normal distribution, we do the same here.  Notice that we     #
               # don't do any weighted mean, and this is because we are looking at the     #
               # species distribution of mortality rates.                                  #
               #---------------------------------------------------------------------------#
               sel.pft           = nplantpft > 0
               if (any(sel.pft)){
                  #------ Save LAI, AGB, and BA for the first time only. ------------------#
                  if (m == muse){
                     ts.lai.plot[m,u] = sum(x=laipft[sel.pft])
                     ts.ba.plot [m,u] = sum(x=bapft [sel.pft])
                     ts.agb.plot[m,u] = sum(x=agbpft[sel.pft])
                  }else{
                     ts.lai.plot[m,u] = ts.lai.plot[muse,u]
                     ts.ba.plot [m,u] = ts.ba.plot [muse,u]
                     ts.agb.plot[m,u] = ts.agb.plot[muse,u]
                  }#end if
                  #------------------------------------------------------------------------#
                  ts.growth.plot[m,u] = sum(x=growthpft[sel.pft]) / sum(nplantpft[sel.pft])
                  ts.mort.plot  [m,u] = weighted.mean( x = mortpft  [sel.pft]
                                                     , w = nplantpft[sel.pft])
                  ts.dimort.plot[m,u] = weighted.mean( x = dimortpft[sel.pft]
                                                     , w = nplantpft[sel.pft])
                  ts.ddmort.plot[m,u] = weighted.mean( x = ddmortpft[sel.pft]
                                                     , w = nplantpft[sel.pft])
                  ts.recr       [m,u] = weighted.mean( x = recrpft  [sel.pft]
                                                     , w = nplantpft[sel.pft])
               }else{
                  ts.growth.plot[m,u] = NA
                  ts.mort.plot  [m,u] = NA
                  ts.dimort.plot[m,u] = NA
                  ts.ddmort.plot[m,u] = NA
                  ts.recr       [m,u] = NA
               }#end if
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Build the size (DBH) structure arrays.                                #
               #---------------------------------------------------------------------------#
               nplantpftdbh = matrix(0.,nrow=n.dbh,ncol=npft)
               agbpftdbh    = matrix(0.,nrow=n.dbh,ncol=npft)
               bapftdbh     = matrix(0.,nrow=n.dbh,ncol=npft)
               laipftdbh    = matrix(0.,nrow=n.dbh,ncol=npft)
               mortpftdbh   = matrix(NA,nrow=n.dbh,ncol=npft)
               dimortpftdbh = matrix(NA,nrow=n.dbh,ncol=npft)
               ddmortpftdbh = matrix(NA,nrow=n.dbh,ncol=npft)
               growthpftdbh = matrix(NA,nrow=n.dbh,ncol=npft)
               for (d in 1:n.dbh){
                  if (all(is.na(dbhfac))){
                     seldbh  = rep(FALSE,times=length(dbhfac))
                  }else{
                     seldbh  = dbhfac == d
                  }#end if
                  for (p in 1:npft){
                     selpft   = pftconow == p
                     sel.dbh  = selpft & seldbh & censtatusconow  > 0
                     sel.cs2  = selpft & seldbh & censtatusconow == 2

                     if (any(sel.dbh)){
                        agbpftdbh    [d,p] = sum( w.nplant    [sel.dbh]
                                                * agbconow    [sel.dbh] )
                        bapftdbh     [d,p] = sum( w.nplant    [sel.dbh]
                                                * baconow     [sel.dbh] )
                        laipftdbh    [d,p] = sum( w.lai       [sel.dbh] )
                     }#end if

                     if (any(sel.cs2)){
                        nplantpftdbh [d,p] = sum( w.nplant [sel.cs2])

                        growthpftdbh [d,p] = sum( w.nplant    [sel.cs2]
                                                * growthconow [sel.cs2] )

                        #---- This is the number of survivors and living before. ----------#
                        survivor           = sum( w.nplant[sel.cs2]         )
                        previous           = sum( w.nplant[sel.cs2] 
                                                * exp(mortconow[sel.cs2])   )
                        mortpftdbh   [d,p] = log( previous / survivor )

                        survivor           = sum( w.nplant[sel.cs2]         )
                        previous           = sum( w.nplant[sel.cs2] 
                                                * exp(dimortconow[sel.cs2]) )
                        dimortpftdbh [d,p] = log( previous / survivor )

                        survivor           = sum( w.nplant[sel.cs2] )
                        previous           = sum( w.nplant[sel.cs2]
                                                * exp(ddmortconow[sel.cs2]) )
                        ddmortpftdbh [d,p] = log( previous / survivor )
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for PFT
                  #------------------------------------------------------------------------#
               }#end for DBH
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Find the mean rates amongst all PFTs.                                #
               #---------------------------------------------------------------------------#
               for (d in 1:n.dbh){
                  sel.pft           = nplantpftdbh[d,] > 0
                  if (any(sel.pft)){
                     if (m == muse){
                        ts.agb.size   [d,m,u] = sum(agbpftdbh[d,sel.pft])
                        ts.ba.size    [d,m,u] = sum(bapftdbh [d,sel.pft])
                        ts.lai.size   [d,m,u] = sum(laipftdbh[d,sel.pft])
                     }else{
                        ts.agb.size   [d,m,u] = ts.agb.size   [d,muse,u]
                        ts.ba.size    [d,m,u] = ts.ba.size    [d,muse,u]
                        ts.lai.size   [d,m,u] = ts.lai.size   [d,muse,u]
                     }#end if
                     ts.growth.size[d,m,u] = ( sum(growthpftdbh[d,sel.pft])
                                             / sum(nplantpftdbh[d,sel.pft]) )
                     ts.mort.size  [d,m,u] = weighted.mean(x = mortpftdbh  [d,sel.pft]
                                                          ,w = nplantpftdbh[d,sel.pft])
                     ts.dimort.size[d,m,u] = weighted.mean(x = dimortpftdbh[d,sel.pft]
                                                          ,w = nplantpftdbh[d,sel.pft])
                     ts.ddmort.size[d,m,u] = weighted.mean(x = ddmortpftdbh[d,sel.pft]
                                                          ,w = nplantpftdbh[d,sel.pft])
                  }else{
                     ts.agb.size   [d,m,u] = 0.
                     ts.ba.size    [d,m,u] = 0.
                     ts.lai.size   [d,m,u] = 0.
                     ts.growth.size[d,m,u] = NA
                     ts.mort.size  [d,m,u] = NA
                     ts.dimort.size[d,m,u] = NA
                     ts.ddmort.size[d,m,u] = NA
                  }#end if
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }else{
               print (paste("     * ",basename(myfile)," wasn't found, skipping it..."
                           ,sep=""))
               
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
         #----- Make the point and legend name sequence. ----------------------------------#
         yeara            = c(qapply(X=ts.census.year,INDEX=census.idx,DIM=1
                                    ,FUN=min,na.rm=TRUE))
         yearz            = c(qapply(X=ts.census.year,INDEX=census.idx,DIM=1
                                    ,FUN=max,na.rm=TRUE))
         pspace$pch       = match(yearz,eft.year)
         pspace$leg.label = paste(sort(unique(yeara)),sort(unique(yearz)),sep="-")
         pspace$leg.pch   = match(sort(unique(yearz)),eft.year)
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
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         pspace$paw              = c(qapply( X     = ts.paw
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         pspace$smpot            = c(qapply( X     = ts.smpot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         #---------------------------------------------------------------------------------#



         #----- Plot-level variables. -----------------------------------------------------#
         pspace$lai.plot.mean    = c(qapply( X     = ts.lai.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean 
                                           , na.rm = TRUE
                                           ))
         pspace$ba.plot.mean     = c(qapply( X     = ts.ba.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean 
                                           , na.rm = TRUE
                                           ))
         pspace$agb.plot.mean    = c(qapply( X     = ts.agb.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean 
                                           , na.rm = TRUE
                                           ))
         pspace$recr.mean        = c(qapply( X     = ts.recr
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean 
                                           , na.rm = TRUE
                                           ))
         pspace$mort.plot.mean   = c(qapply( X     = ts.mort.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         pspace$dimort.plot.mean = c(qapply( X     = ts.dimort.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         pspace$ddmort.plot.mean = c(qapply( X     = ts.ddmort.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         pspace$growth.plot.mean = c(qapply( X     = ts.growth.plot
                                           , INDEX = census.idx
                                           , DIM   = 1
                                           , FUN   = mean
                                           , na.rm = TRUE
                                           ))
         #---------------------------------------------------------------------------------#



         #----- Size-level variables. -----------------------------------------------------#
         pspace$lai.size.mean    = matrix(qapply( X     = ts.lai.size
                                                , INDEX = census.idx
                                                , DIM   = 2
                                                , FUN   = mean
                                                , na.rm = TRUE
                                                )
                                         , nrow = n.dbh
                                         , ncol = nfullcyc*(n.census-1)
                                         )#end matrix
         pspace$ba.size.mean     = matrix(qapply( X     = ts.ba.size
                                                , INDEX = census.idx
                                                , DIM   = 2
                                                , FUN   = mean
                                                , na.rm = TRUE
                                                )
                                         , nrow = n.dbh
                                         , ncol = nfullcyc*(n.census-1)
                                         )#end matrix
         pspace$agb.size.mean    = matrix(qapply( X     = ts.agb.size
                                                , INDEX = census.idx
                                                , DIM   = 2
                                                , FUN   = mean
                                                , na.rm = TRUE
                                                )
                                         , nrow = n.dbh
                                         , ncol = nfullcyc*(n.census-1)
                                         )#end matrix
         pspace$mort.size.mean   = matrix(qapply( X     = ts.mort.size
                                                , INDEX = census.idx
                                                , DIM   = 2
                                                , FUN   = mean
                                                , na.rm = TRUE
                                                )
                                         , nrow = n.dbh
                                         , ncol = nfullcyc*(n.census-1)
                                         )#end matrix
         pspace$dimort.size.mean = matrix(qapply( X     = ts.dimort.size
                                                , INDEX = census.idx
                                                , DIM   = 2
                                                , FUN   = mean
                                                , na.rm = TRUE
                                                )
                                         , nrow = n.dbh
                                         , ncol = nfullcyc*(n.census-1)
                                         )#end matrix
         pspace$ddmort.size.mean = matrix(qapply( X     = ts.ddmort.size
                                                , INDEX = census.idx
                                                , DIM   = 2
                                                , FUN   = mean
                                                , na.rm = TRUE
                                                )
                                         , nrow = n.dbh
                                         , ncol = nfullcyc*(n.census-1)
                                         )#end matrix
         pspace$growth.size.mean = matrix(qapply( X     = ts.growth.size
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
      #     Remove the 3rd. dimension for the means.                                       #
      #------------------------------------------------------------------------------------#
      ms.recr.mean        = apply(X=ts.recr       ,MARGIN=1     ,FUN=mean    ,na.rm=TRUE)
      ms.mort.plot.mean   = apply(X=ts.mort.plot  ,MARGIN=1     ,FUN=mean    ,na.rm=TRUE)
      ms.dimort.plot.mean = apply(X=ts.dimort.plot,MARGIN=1     ,FUN=mean    ,na.rm=TRUE)
      ms.ddmort.plot.mean = apply(X=ts.ddmort.plot,MARGIN=1     ,FUN=mean    ,na.rm=TRUE)
      ms.growth.plot.mean = apply(X=ts.growth.plot,MARGIN=1     ,FUN=mean    ,na.rm=TRUE)
      ms.mort.size.mean   = apply(X=ts.mort.size  ,MARGIN=c(1,2),FUN=mean    ,na.rm=TRUE)
      ms.dimort.size.mean = apply(X=ts.dimort.size,MARGIN=c(1,2),FUN=mean    ,na.rm=TRUE)
      ms.ddmort.size.mean = apply(X=ts.ddmort.size,MARGIN=c(1,2),FUN=mean    ,na.rm=TRUE)
      ms.growth.size.mean = apply(X=ts.growth.size,MARGIN=c(1,2),FUN=mean    ,na.rm=TRUE)
      ms.recr.mean        [! is.finite(ms.recr.mean       )] = NA
      ms.mort.plot.mean   [! is.finite(ms.mort.plot.mean  )] = NA
      ms.dimort.plot.mean [! is.finite(ms.dimort.plot.mean)] = NA
      ms.ddmort.plot.mean [! is.finite(ms.ddmort.plot.mean)] = NA
      ms.growth.plot.mean [! is.finite(ms.growth.plot.mean)] = NA
      ms.mort.size.mean   [! is.finite(ms.mort.size.mean  )] = NA
      ms.dimort.size.mean [! is.finite(ms.dimort.size.mean)] = NA
      ms.ddmort.size.mean [! is.finite(ms.ddmort.size.mean)] = NA
      ms.growth.size.mean [! is.finite(ms.growth.size.mean)] = NA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the average rates for the census period using log-normal.                 #
      #------------------------------------------------------------------------------------#
      print(paste("   - Finding the average rates...",sep=""))
      model$cwhen     = census.obs$when
      #----- Plot-level. ------------------------------------------------------------------#
      model$recr               = list()
      model$recr$mean          = c(NA,tapply( X     = ms.recr.mean
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$recr$median        = c(NA,tapply( X     = ms.recr.mean
                                            , INDEX = census.idx
                                            , FUN   = median
                                            , na.rm = TRUE
                                            ) )
      model$recr$q025          = c(NA,tapply( X     = ms.recr.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.025
                                            , na.rm = TRUE
                                            ) )
      model$recr$q975          = c(NA,tapply( X     = ms.recr.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.975
                                            , na.rm = TRUE
                                            ) )
      model$mort.plot          = list()
      model$mort.plot$mean     = c(NA,tapply( X     = ms.mort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$mort.plot$median   = c(NA,tapply( X     = ms.mort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = median
                                            , na.rm = TRUE
                                            ) )
      model$mort.plot$q025     = c(NA,tapply( X     = ms.mort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.025
                                            , na.rm = TRUE
                                            ) )
      model$mort.plot$q975     = c(NA,tapply( X     = ms.mort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.975
                                            , na.rm = TRUE
                                            ) )
      model$ddmort.plot        = list()
      model$ddmort.plot$mean   = c(NA,tapply( X     = ms.ddmort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$ddmort.plot$median = c(NA,tapply( X     = ms.ddmort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = median
                                            , na.rm = TRUE
                                            ) )
      model$ddmort.plot$q025   = c(NA,tapply( X     = ms.ddmort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.025
                                            , na.rm = TRUE
                                            ) )
      model$ddmort.plot$q975   = c(NA,tapply( X     = ms.ddmort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.975
                                            , na.rm = TRUE
                                            ) )
      model$dimort.plot        = list()
      model$dimort.plot$mean   = c(NA,tapply( X     = ms.dimort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$dimort.plot$median = c(NA,tapply( X     = ms.dimort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = median
                                            , na.rm = TRUE
                                            ) )
      model$dimort.plot$q025   = c(NA,tapply( X     = ms.dimort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.025
                                            , na.rm = TRUE
                                            ) )
      model$dimort.plot$q975   = c(NA,tapply( X     = ms.dimort.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.975
                                            , na.rm = TRUE
                                            ) )
      model$growth.plot        = list()
      model$growth.plot$mean   = c(NA,tapply( X     = ms.growth.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$growth.plot$median = c(NA,tapply( X     = ms.growth.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = median
                                            , na.rm = TRUE
                                            ) )
      model$growth.plot$q025   = c(NA,tapply( X     = ms.growth.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.025
                                            , na.rm = TRUE
                                            ) )
      model$growth.plot$q975   = c(NA,tapply( X     = ms.growth.plot.mean
                                            , INDEX = census.idx
                                            , FUN   = quantile
                                            , probs = 0.975
                                            , na.rm = TRUE
                                            ) )
      #----- Size- and plot-level. --------------------------------------------------------#
      empty           = rep(NA,times=n.dbh)
      model$mort.size          = list()
      model$mort.size$mean     = cbind(empty,qapply( X     = ms.mort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = mean
                                                   , na.rm = TRUE
                                                   ) )
      model$mort.size$median   = cbind(empty,qapply( X     = ms.mort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = median
                                                   , na.rm = TRUE
                                                   ) )
      model$mort.size$q025     = cbind(empty,qapply( X     = ms.mort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.025
                                                   , na.rm = TRUE
                                                   ) )
      model$mort.size$q975     = cbind(empty,qapply( X     = ms.mort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.975
                                                   , na.rm = TRUE
                                                   ) )
      model$ddmort.size        = list()
      model$ddmort.size$mean   = cbind(empty,qapply( X     = ms.ddmort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = mean
                                                   , na.rm = TRUE
                                                   ) )
      model$ddmort.size$median = cbind(empty,qapply( X     = ms.ddmort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = median
                                                   , na.rm = TRUE
                                                   ) )
      model$ddmort.size$q025   = cbind(empty,qapply( X     = ms.ddmort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.025
                                                   , na.rm = TRUE
                                                   ) )
      model$ddmort.size$q975   = cbind(empty,qapply( X     = ms.ddmort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.975
                                                   , na.rm = TRUE
                                                   ) )
      model$dimort.size        = list()
      model$dimort.size$mean   = cbind(empty,qapply( X     = ms.dimort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = mean
                                                   , na.rm = TRUE
                                                   ) )
      model$dimort.size$median = cbind(empty,qapply( X     = ms.dimort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = median
                                                   , na.rm = TRUE
                                                   ) )
      model$dimort.size$q025   = cbind(empty,qapply( X     = ms.dimort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.025
                                                   , na.rm = TRUE
                                                   ) )
      model$dimort.size$q975   = cbind(empty,qapply( X     = ms.dimort.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.975
                                                   , na.rm = TRUE
                                                   ) )
      model$growth.size        = list()
      model$growth.size$mean   = cbind(empty,qapply( X     = ms.growth.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = mean
                                                   , na.rm = TRUE
                                                   ) )
      model$growth.size$median = cbind(empty,qapply( X     = ms.growth.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = median
                                                   , na.rm = TRUE
                                                   ) )
      model$growth.size$q025   = cbind(empty,qapply( X     = ms.growth.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.025
                                                   , na.rm = TRUE
                                                   ) )
      model$growth.size$q975   = cbind(empty,qapply( X     = ms.growth.size.mean
                                                   , DIM   = 2
                                                   , INDEX = census.idx
                                                   , FUN   = quantile
                                                   , probs = 0.975
                                                   , na.rm = TRUE
                                                   ) )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     In case vegetation dynamics is turned off, all variables above will be NA.     #
      #------------------------------------------------------------------------------------#
      model$recr$mean          [!is.finite(model$recr$mean         )] = 0.
      model$recr$median        [!is.finite(model$recr$median       )] = 0.
      model$recr$q025          [!is.finite(model$recr$q025         )] = 0.
      model$recr$q975          [!is.finite(model$recr$q975         )] = 0.
      model$mort.plot$mean     [!is.finite(model$mort.plot$mean    )] = 0.
      model$mort.plot$median   [!is.finite(model$mort.plot$median  )] = 0.
      model$mort.plot$q025     [!is.finite(model$mort.plot$q025    )] = 0.
      model$mort.plot$q975     [!is.finite(model$mort.plot$q975    )] = 0.
      model$ddmort.plot$mean   [!is.finite(model$ddmort.plot$mean  )] = 0.
      model$ddmort.plot$median [!is.finite(model$ddmort.plot$median)] = 0.
      model$ddmort.plot$q025   [!is.finite(model$ddmort.plot$q025  )] = 0.
      model$ddmort.plot$q975   [!is.finite(model$ddmort.plot$q975  )] = 0.
      model$dimort.plot$mean   [!is.finite(model$dimort.plot$mean  )] = 0.
      model$dimort.plot$median [!is.finite(model$dimort.plot$median)] = 0.
      model$dimort.plot$q025   [!is.finite(model$dimort.plot$q025  )] = 0.
      model$dimort.plot$q975   [!is.finite(model$dimort.plot$q975  )] = 0.
      model$growth.plot$mean   [!is.finite(model$growth.plot$mean  )] = 0.
      model$growth.plot$median [!is.finite(model$growth.plot$median)] = 0.
      model$growth.plot$q025   [!is.finite(model$growth.plot$q025  )] = 0.
      model$growth.plot$q975   [!is.finite(model$growth.plot$q975  )] = 0.
      #----- Size- and plot-level. --------------------------------------------------------#
      model$mort.size$mean     [!is.finite(model$mort.size$mean    )] = 0.
      model$mort.size$median   [!is.finite(model$mort.size$median  )] = 0.
      model$mort.size$q025     [!is.finite(model$mort.size$q025    )] = 0.
      model$mort.size$q975     [!is.finite(model$mort.size$q975    )] = 0.
      model$ddmort.size$mean   [!is.finite(model$ddmort.size$mean  )] = 0.
      model$ddmort.size$median [!is.finite(model$ddmort.size$median)] = 0.
      model$ddmort.size$q025   [!is.finite(model$ddmort.size$q025  )] = 0.
      model$ddmort.size$q975   [!is.finite(model$ddmort.size$q975  )] = 0.
      model$dimort.size$mean   [!is.finite(model$dimort.size$mean  )] = 0.
      model$dimort.size$median [!is.finite(model$dimort.size$median)] = 0.
      model$dimort.size$q025   [!is.finite(model$dimort.size$q025  )] = 0.
      model$dimort.size$q975   [!is.finite(model$dimort.size$q975  )] = 0.
      model$growth.size$mean   [!is.finite(model$growth.size$mean  )] = 0.
      model$growth.size$median [!is.finite(model$growth.size$median)] = 0.
      model$growth.size$q025   [!is.finite(model$growth.size$q025  )] = 0.
      model$growth.size$q975   [!is.finite(model$growth.size$q975  )] = 0.
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Make the directories.                                                          #
      #------------------------------------------------------------------------------------#
      outplot   = paste(outpref,"census_plot"  ,sep="/")
      outsize   = paste(outpref,"census_size"  ,sep="/")
      if (! file.exists(outplot)) dir.create(outplot)
      if (! file.exists(outsize)) dir.create(outsize)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the time series of all variables.                                        #
      #------------------------------------------------------------------------------------#
      print(paste("   - Plotting the plot-level rates...",sep=""))
      for (p in 1:nplotvar){
         #---------------------------------------------------------------------------------#
         #      Copy to scratch variables.                                                 #
         #---------------------------------------------------------------------------------#
         this.plot       = plotvar[[p]]
         this.vnam.ed    = this.plot$vnam.ed
         this.vnam.obs   = this.plot$vnam.obs
         this.desc       = this.plot$desc
         this.unit       = this.plot$unit
         this.col.obser  = this.plot$col.obser
         this.col.model  = this.plot$col.model
         this.leg.corner = this.plot$leg.corner
         this.plog       = this.plot$plog
         if (this.plog){
            plog = "y"
         }else{
            plog =""
         }#end if
         print(paste("     * ",this.desc,"...",sep=""))
         #---------------------------------------------------------------------------------#

         #----- Load the data. ------------------------------------------------------------#
         when            = census.obs$when[2:n.census]
         this.obs.mean   = 100. * census.obs[[this.vnam.obs]]       [4,2:n.census]
         this.obs.q025   = 100. * census.obs[[this.vnam.obs]]       [5,2:n.census]
         this.obs.q975   = 100. * census.obs[[this.vnam.obs]]       [6,2:n.census]
         this.mod.mean   = 100. * model     [[this.vnam.ed ]]$median[  2:n.census]
         this.mod.q025   = 100. * model     [[this.vnam.ed ]]$q025  [  2:n.census]
         this.mod.q975   = 100. * model     [[this.vnam.ed ]]$q975  [  2:n.census]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the x axis.                                                           #
         #---------------------------------------------------------------------------------#
         whenplot = pretty.time(when,n=6)
         #---------------------------------------------------------------------------------#

         #----- Find the plot range. ------------------------------------------------------#
         ylim.test = c(this.obs.q025,this.obs.q975,this.mod.q025,this.mod.q975)
         ylim.test = ylim.test[is.finite(ylim.test) & ylim.test > 0 ]
         ylimit    = range(ylim.test,na.rm=TRUE)
         #----- Make room for the legend. -------------------------------------------------#
         if (this.plog) ylimit = log(ylimit)
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
         if (this.plog) ylimit = exp(ylimit)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make the error polygon.                                                     #
         #---------------------------------------------------------------------------------#
         err.x = c(when,rev(when),NA,when,rev(when))
         err.y = c(this.obs.q025,rev(this.obs.q975),NA,this.mod.q025,rev(this.mod.q975))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            #----- Open file. -------------------------------------------------------------#
            fichier = paste(outplot,"/",this.vnam.ed,"-",suffix,".",outform[o],sep="")
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
            #------------------------------------------------------------------------------#


            #----- Make plot annotation. --------------------------------------------------#
            letitre = paste(lieu,"\n","Average ",this.desc,sep="")
            lex     = paste("Census time")
            ley     = paste(this.desc,this.unit,sep=" ")
            #------------------------------------------------------------------------------#


            #----- Start the plot. --------------------------------------------------------#
            plot(x=when,y=this.obs.mean,type="n",main=letitre,xlab=lex,ylab=ley,ylim=ylimit
                ,log=plog,xaxt="n",cex.main=cex.main)
            #----- Special, time-friendly X-Axis and grid. --------------------------------#
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){
               abline(v=whenplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            #----- Confidence interval. ---------------------------------------------------#
            epolygon(x=err.x,y=err.y,col=c(this.col.obser[2],this.col.model[2])
                    ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
            #----- Observed data. ---------------------------------------------------------#
            points(x=when,y=this.obs.mean,col=this.col.obser[1],lwd=3.0
                  ,type="o",pch=16,cex=1.0)
            #----- Modelled data. ---------------------------------------------------------#
            points(x=when,y=this.mod.mean,col=this.col.model[1],lwd=3.0
                  ,type="o",pch=16,cex=1.0)
            #----- Legend. ----------------------------------------------------------------#
            legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                  ,fill   = c(this.col.obser[2],this.col.model[2])
                  ,border = c(this.col.obser[2],this.col.model[2])
                  ,angle=c(-45,45),density=60,lwd=3.0
                  ,col=c(this.col.obser[1],this.col.model[1])
                  ,bg="white",title="Shaded area - 95%C.I.",cex=1.0,pch=16)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Close the device.                                                        #
            #------------------------------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#
         }#end for (o in 1:nout)
      }#end for (p in 1:nplotvar)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the time series of all variables.                                        #
      #------------------------------------------------------------------------------------#
      print(paste("   - Plotting the DBH-dependent plot-level rates...",sep=""))
      for (p in 1:nsizevar){
         #---------------------------------------------------------------------------------#
         #      Copy to scratch variables.                                                 #
         #---------------------------------------------------------------------------------#
         this.size       = sizevar[[p]]
         this.vnam.ed    = this.size$vnam.ed
         this.vnam.obs   = this.size$vnam.obs
         this.desc       = this.size$desc
         this.unit       = this.size$unit
         this.col.obser  = this.size$col.obser
         this.col.model  = this.size$col.model
         this.leg.corner = this.size$leg.corner
         this.plog       = this.size$plog
         if (this.plog){
            plog = "y"
         }else{
            plog =""
         }#end if
         print(paste("     * ",this.desc,"...",sep=""))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Create one directory per variable.                                         #
         #---------------------------------------------------------------------------------#
         outvar = paste(outsize,this.vnam.ed,sep="/")
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#



         #----- Load the data. ------------------------------------------------------------#
         when            = census.obs$when[2:n.census]
         this.obs.mean   = 100. * census.obs[[this.vnam.obs]]       [4,,2:n.census]
         this.obs.q025   = 100. * census.obs[[this.vnam.obs]]       [5,,2:n.census]
         this.obs.q975   = 100. * census.obs[[this.vnam.obs]]       [6,,2:n.census]
         this.mod.mean   = 100. * model     [[this.vnam.ed ]]$median[  ,2:n.census]
         this.mod.q025   = 100. * model     [[this.vnam.ed ]]$q025  [  ,2:n.census]
         this.mod.q975   = 100. * model     [[this.vnam.ed ]]$q975  [  ,2:n.census]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Find the x axis.                                                           #
         #---------------------------------------------------------------------------------#
         whenplot = pretty.time(when,n=6)
         #---------------------------------------------------------------------------------#

         #----- Find the plot range. ------------------------------------------------------#
         ylim.test = c(this.obs.q025,this.obs.q975,this.mod.q025,this.mod.q975)
         ylim.test = ylim.test[is.finite(ylim.test) & ylim.test > 0 ]
         ylimit    = range(ylim.test,na.rm=TRUE)
         #----- Make room for the legend. -------------------------------------------------#
         if (this.plog) ylimit = log(ylimit)
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
         if (this.plog) ylimit = exp(ylimit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop over DBH classes.                                                      #
         #---------------------------------------------------------------------------------#
         for (d in 1:n.dbh){
            print(paste("       ~ ",dbh.names[d],"...",sep=""))

            #------------------------------------------------------------------------------#
            #     Make the error polygon.                                                  #
            #------------------------------------------------------------------------------#
            err.x = c(when,rev(when),NA,when,rev(when))
            err.y = c(this.obs.q025[d,],rev(this.obs.q975[d,])
                     ,NA
                     ,this.mod.q025[d,],rev(this.mod.q975[d,]))
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Make the label.                                                          #
            #------------------------------------------------------------------------------#
            if (d == 1){
               dbh.low   = sprintf("%3.3i",0)
            }else{
               dbh.low   = sprintf("%3.3i",census.obs$dbh.breaks[  d])
            }#end if
            if (d == n.dbh){
               dbh.high  = "Inf"
            }else{
               dbh.high  = sprintf("%3.3i",census.obs$dbh.breaks[d+1])
            }#end if
            dbh.label = paste(dbh.low,dbh.high,sep="-")
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open file. ----------------------------------------------------------#
               fichier = paste(outvar,"/",this.vnam.ed,"-dbh_",dbh.label
                              ,"-",suffix,".",outform[o],sep="")
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


               #----- Make plot annotation. -----------------------------------------------#
               letitre = paste(lieu,"\n","Average ",this.desc
                                        ," - DBH Class: ",dbh.names[d],sep="")
               lex     = paste("Census time")
               ley     = paste(this.desc,this.unit,sep=" ")
               #---------------------------------------------------------------------------#


               #----- Start the plot. -----------------------------------------------------#
               plot(x=when,y=this.obs.mean[d,],type="n",main=letitre,xlab=lex,ylab=ley
                   ,ylim=ylimit,log=plog,xaxt="n",cex.main=cex.main)
               #----- Special, time-friendly X-Axis and grid. -----------------------------#
               axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
               if (plotgrid){
                  abline(v=whenplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
               }#end if
               #----- Confidence interval. ------------------------------------------------#
               epolygon(x=err.x,y=err.y,col=c(this.col.obser[2],this.col.model[2])
                       ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
               #----- Observed data. ------------------------------------------------------#
               points(x=when,y=this.obs.mean[d,],col=this.col.obser[1],lwd=3.0
                     ,type="o",pch=16,cex=1.0)
               #----- Modelled data. ------------------------------------------------------#
               points(x=when,y=this.mod.mean[d,],col=this.col.model[1],lwd=3.0
                     ,type="o",pch=16,cex=1.0)
               #----- Legend. -------------------------------------------------------------#
               legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                     ,fill   = c(this.col.obser[2],this.col.model[2])
                     ,border = c(this.col.obser[2],this.col.model[2])
                     ,angle=c(-45,45),density=40,lwd=3.0
                     ,col=c(this.col.obser[1],this.col.model[1])
                     ,bg="white",title="Shaded area - 95%C.I.",cex=1.0,pch=16)
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
            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end for (d in 1:n.dbh)
         #---------------------------------------------------------------------------------#
      }#end for (p in 1:nplotvar)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Plot the time series of all variables.                                        #
      #------------------------------------------------------------------------------------#
      print(paste("   - Plotting the structure of DBH-dependent rates...",sep=""))
      for (p in 1:nsizevar){
         #---------------------------------------------------------------------------------#
         #      Copy to scratch variables.                                                 #
         #---------------------------------------------------------------------------------#
         this.size       = sizevar[[p]]
         this.vnam.ed    = this.size$vnam.ed
         this.vnam.obs   = this.size$vnam.obs
         this.desc       = this.size$desc
         this.unit       = this.size$unit
         this.col.obser  = this.size$col.obser
         this.col.model  = this.size$col.model
         this.leg.corner = this.size$leg.corner
         this.plog       = this.size$plog
         if (this.plog){
            plog = "y"
         }else{
            plog =""
         }#end if
         print(paste("     * ",this.desc,"...",sep=""))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Create one directory per variable.                                         #
         #---------------------------------------------------------------------------------#
         outvar = paste(outsize,this.vnam.ed,sep="/")
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#



         #----- Load the data. ------------------------------------------------------------#
         yr.census       = numyears(census.obs$when)
         this.obs.mean  = 100. * census.obs[[this.vnam.obs]]       [4,,] 
         this.obs.q025  = 100. * census.obs[[this.vnam.obs]]       [5,,] 
         this.obs.q975  = 100. * census.obs[[this.vnam.obs]]       [6,,] 
         this.mod.mean  = 100. * model     [[this.vnam.ed ]]$median[  ,]
         this.mod.q025  = 100. * model     [[this.vnam.ed ]]$q025  [  ,]
         this.mod.q975  = 100. * model     [[this.vnam.ed ]]$q975  [  ,]
         #---------------------------------------------------------------------------------#


         #----- Find the plot range. ------------------------------------------------------#
         ylim.test = c(this.obs.q025,this.obs.q975,this.mod.q025,this.mod.q975)
         ylim.test = ylim.test[is.finite(ylim.test) & ylim.test > 0 ]
         ylimit    = range(ylim.test,na.rm=TRUE)
         #----- Make room for the legend. -------------------------------------------------#
         if (this.plog) ylimit = log(ylimit)
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
         if (this.plog) ylimit = exp(ylimit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop over DBH classes.                                                      #
         #---------------------------------------------------------------------------------#
         for (y in 2:n.census){
            period = paste("(",yr.census[y-1],"-",yr.census[y],")",sep="")
            print(paste("       ~ Period: ",period,"...",sep=""))

            #------------------------------------------------------------------------------#
            #     Make the error polygon.                                                  #
            #------------------------------------------------------------------------------#
            err.x = c(x.dbh,rev(x.dbh),NA,x.dbh,rev(x.dbh))
            err.y = c(this.obs.q025[,y],rev(this.obs.q975[,y])
                     ,NA
                     ,this.mod.q025[,y],rev(this.mod.q975[,y]))
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Make the label.                                                          #
            #------------------------------------------------------------------------------#
            yr.label = paste(yr.census[y-1],yr.census[y],sep="-")
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Loop over formats.                                                      #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               #----- Open file. ----------------------------------------------------------#
               fichier = paste(outvar,"/",this.vnam.ed,"-census",yr.label
                              ,"-",suffix,".",outform[o],sep="")
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


               #----- Make plot annotation. -----------------------------------------------#
               letitre = paste(lieu,"\n","Average ",this.desc," - Period: ",period,sep="")
               lex     = paste("Minimum DBH for this class [cm]")
               ley     = paste(this.desc,this.unit,sep=" ")
               #---------------------------------------------------------------------------#


               #----- Start the plot. -----------------------------------------------------#
               plot(x=x.dbh,y=this.obs.mean[,y],type="n",main=letitre,xlab=lex,ylab=ley
                   ,ylim=ylimit,log=plog,cex.main=cex.main)
               #----- Plot the grid. ------------------------------------------------------#
               if (plotgrid) grid(col="gray52",lty="solid")
               #----- Confidence interval. ------------------------------------------------#
               epolygon(x=err.x,y=err.y,col=c(this.col.obser[2],this.col.model[2])
                       ,angle=c(-45,45),density=40,lty="solid",lwd=1.0)
               #----- Observed data. ------------------------------------------------------#
               points(x=x.dbh,y=this.obs.mean[,y],col=this.col.obser[1],lwd=3.0,type="o"
                     ,pch=16,cex=1.0)
               #----- Modelled data. ------------------------------------------------------#
               points(x=x.dbh,y=this.mod.mean[,y],col=this.col.model[1],lwd=3.0,type="o"
                     ,pch=16,cex=1.0)
               #----- Legend. -------------------------------------------------------------#
               legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                     ,fill   = c(this.col.obser[2],this.col.model[2])
                     ,border = c(this.col.obser[2],this.col.model[2])
                     ,angle=c(-45,45),density=40,lwd=3.0
                     ,col=c(this.col.obser[1],this.col.model[1])
                     ,bg="white",title="Shaded area - 95%C.I.",cex=1.0,pch=16)
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
            }#end for (o in 1:nout)
            #------------------------------------------------------------------------------#
         }#end for (d in 1:n.dbh)
         #---------------------------------------------------------------------------------#
      }#end for (p in 1:nplotvar)
      #------------------------------------------------------------------------------------#
   }#end if (census.name %in% ls())
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #=======================================================================================#
   #     Plot the monthly mean variables as functions of other 2 environment variables.    #
   #---------------------------------------------------------------------------------------#
   print(paste(" + Plotting parameter space...",sep=""))


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
      print(paste("     * Y: ",ydesc,"..."))
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

      sel                   = is.finite(yvar.plot) & ( yvar.plot > 0 | (! ylog))
      ylimit.plot           = range(yvar.plot[sel],na.rm=TRUE)
      if (ylog) ylimit.plot = log(ylimit.plot)
      ylimit.plot[2]        = ylimit.plot[2] + scalleg * diff(ylimit.plot)
      if (ylog) ylimit.plot = exp(ylimit.plot)

      if (sizetoo){
         ylimit.size = matrix(nrow=n.dbh,ncol=2)
         for (d in 1:n.dbh){
            sel                       = ( is.finite(yvar.size[d,]) 
                                        & ( yvar.size[d,] > 0 | (! ylog)) )
            ylimit.size[d,]           = range(yvar.size[d,sel],na.rm=TRUE)
            if (ylog) ylimit.size[d,] = log(ylimit.size[d,])
            ylimit.size[d,2]          = ( ylimit.size[d,2]
                                        + scalleg * diff(ylimit.size[d,]) )
            if (ylog) ylimit.size[d,] = exp(ylimit.size[d,])
         }#end for
      }else{
         ylimit.size = ylimit.plot
      }#end if
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
            zlog     = this.z$log
            print(paste("      ~ X:",xdesc,"   Z: ",zdesc,"..."))
            #------------------------------------------------------------------------------#



            #----- Annotation for the colour map ("Z" axis). ------------------------------#
            zvar  = zmult * ( pspace[[zvname]] + zadd )
            lez   = paste(zkey,"\n [",zunit,"]",sep="")
            #------------------------------------------------------------------------------#



            #----- Find the range for the scale. ------------------------------------------#
            sel     = is.finite(zvar) & ( zvar > 0 | (! zlog) )
            zlimit  = range(zvar[sel],na.rm=TRUE)
            #------------------------------------------------------------------------------#


            #----- Title. -----------------------------------------------------------------#
            letitre.plot = paste(lieu,paste(zdesc,"Plot level",sep=" - "),sep="\n")
            letitre.size = paste(lieu,paste(zdesc," - DBH Class:",dbh.names,sep="")
                                ,sep="\n")
            #------------------------------------------------------------------------------#



            #----- Attribute symbols according to the year. -------------------------------#
            this.pch  = pspace$pch
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Make the list of commands for plot.axes.                                 #
            #------------------------------------------------------------------------------#
            plot.axes      = list()
            plot.axes[[1]] = list( x.axis = list(side=1)
                                 , y.axis = list(side=2)
                                 , grid   = list(col="gray62",lty="solid")
                                 , legend = list( x      = leg.pos
                                                , inset  = 0.01
                                                , legend = pspace$leg.label
                                                , col    = "black"
                                                , bg     = "white"
                                                , pch    = pspace$leg.pch
                                                , title  = "Census"
                                                , ncol   = 2
                                                , pt.cex = 1./0.9
                                                , cex    = 0.9
                                                )#end legend
                                 )#end list
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the plot-level.                                                     #
            #------------------------------------------------------------------------------#
            print(paste("         > Plot-level..."))
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
               colourmap(x=xvar.plot,y=yvar.plot,z=zvar
                        ,xlim=xlimit.plot,ylim=ylimit.plot,zlim=zlimit
                        ,colour.palette=muitas,cex=1.6,pch=this.pch,lwd=3,log=plog
                        ,plot.title=title(main=letitre.plot,xlab=lex,ylab=ley
                                         ,cex.main=cex.main)
                        ,key.title=title(main=lez,cex.main=0.8),key.log=zlog
                        ,plot.axes=plot.axes
                        )#end colourmap
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

               #---------------------------------------------------------------------------#
               #    Loop over DBH classes.                                                 #
               #---------------------------------------------------------------------------#
               for (d in 1:n.dbh){
                  print(paste("         > DBH class: ",dbh.names[d],"..."))

                  #------------------------------------------------------------------------#
                  #     Make the label.                                                    #
                  #------------------------------------------------------------------------#
                  if (d == 1){
                     dbh.low   = sprintf("%3.3i",0)
                  }else{
                     dbh.low   = sprintf("%3.3i",census.obs$dbh.breaks[  d])
                  }#end if
                  if (d == n.dbh){
                     dbh.high  = "Inf"
                  }else{
                     dbh.high  = sprintf("%3.3i",census.obs$dbh.breaks[d+1])
                  }#end if
                  dbh.label = paste(dbh.low,dbh.high,sep="-")
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Plot the bar plot.                                                 #
                  #------------------------------------------------------------------------#
                  for (o in 1:nout){
                     #----- Open the file. ------------------------------------------------#
                     fichier = paste(outyvar,"/size",dbh.label,"_x_",xvname
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
                     #---------------------------------------------------------------------#



                     #----- Plot the parameter space. -------------------------------------#
                     colourmap(x=xvar.size[d,],y=yvar.size[d,],z=zvar
                              ,xlim=xlimit.size[d,],ylim=ylimit.size[d,],zlim=zlimit
                              ,colour.palette=muitas,cex=1.6,pch=this.pch,lwd=3,log=plog
                              ,plot.title=title(main=letitre.size[d],cex.main=cex.main
                                               ,xlab=lex,ylab=ley)
                              ,key.title=title(main=lez,cex.main=0.8),key.log=zlog
                              ,plot.axes=plot.axes
                              )#end colourmap
                     #---------------------------------------------------------------------#


                     #---------------------------------------------------------------------#
                     #     Close the device.                                               #
                     #---------------------------------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     #---------------------------------------------------------------------#
                  }#end for outform
                  #------------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#

         }#end for
         #---------------------------------------------------------------------------------#



      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#

#q("no")
