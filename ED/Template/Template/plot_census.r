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
myplaces       = c("thispoly")
sasmonth.short = c(2,5,8,11)
sasmonth.long  = 5
nyears.long    = 25
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

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.20
cex.main       = 0.8            # Scale coefficient for the title
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
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "chartreuse4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 2]] = list( vnam.ed    = "mort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Total mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "purple4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 3]] = list( vnam.ed    = "ddmort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Density dependent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "darkorchid"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 4]] = list( vnam.ed    = "dimort.plot"
                    , vnam.obs   = "mort.plot"
                    , desc       = "Density independent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "purple"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
plotvar[[ 5]] = list( vnam.ed    = "growth.plot"
                    , vnam.obs   = "growth.plot"
                    , desc       = "Growth rate"
                    , unit       = "[%DBH/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "royalblue4"
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
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "purple4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 2]] = list( vnam.ed    = "ddmort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Density dependent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "darkorchid"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 3]] = list( vnam.ed    = "dimort.size"
                    , vnam.obs   = "mort.size"
                    , desc       = "Density independent mortality rate"
                    , unit       = "[%pop/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "purple"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
sizevar[[ 4]] = list( vnam.ed    = "growth.size"
                    , vnam.obs   = "growth.size"
                    , desc       = "Growth rate"
                    , unit       = "[%DBH/yr]"
                    , col.obser  = c("gray42","gray21")
                    , col.model  = "royalblue4"
                    , leg.corner = "topleft"
                    , plog       = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#


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
source(paste(srcdir,"cloudy.r"          ,sep="/"))
source(paste(srcdir,"epolygon.r"       ,sep="/"))
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
      n.months = length(census.idx)
      #------------------------------------------------------------------------------------#



      #----- Initialise all the structures for which we will compare. ---------------------#
      ts.recr          = rep(NA,n.months)
      ts.mort.plot     = rep(NA,n.months)
      ts.ddmort.plot   = rep(NA,n.months)
      ts.dimort.plot   = rep(NA,n.months)
      ts.growth.plot   = rep(NA,n.months)
      ts.mort.size     = matrix(nrow=n.dbh,ncol=n.months)
      ts.ddmort.size   = matrix(nrow=n.dbh,ncol=n.months)
      ts.dimort.size   = matrix(nrow=n.dbh,ncol=n.months)
      ts.growth.size   = matrix(nrow=n.dbh,ncol=n.months)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop over all times, and retrieve the data.                                    #
      #------------------------------------------------------------------------------------#
      now.month = nummonths(census.obs$when[1])
      now.year  = numyears (census.obs$when[1])
      for (m in 1:n.months){
         now.month = (now.month %% 12) + 1
         now.year  = now.year + as.integer(now.month == 1)

         #----- Build the file name. ------------------------------------------------------#
         cmonth    = sprintf("%2.2i",now.month)
         cyear     = sprintf("%2.2i",now.year )
         myfile = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5",sep="")
         #---------------------------------------------------------------------------------#



         #----- Read data and close connection immediately after. -------------------------#
         print (paste("     * Reading ",basename(myfile),"...",sep=""))
         mymont = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
         #---------------------------------------------------------------------------------#


         #---- Read in the site-level area. -----------------------------------------------#
         areasi     = mymont$AREA.SI
         npatches   = mymont$SIPA.N
         #---------------------------------------------------------------------------------#


         #----- Read a few patch-level variables. -----------------------------------------#
         areapa     = mymont$AREA * rep(areasi,times=npatches)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Read the cohort-level variables.  Because empty patchs do exist (deserts),  #
         # we must check whether there is any cohort to be read.  If not, assign NA to all #
         # variables.                                                                      #
         #---------------------------------------------------------------------------------#
         ncohorts   = mymont$PACO.N
         if (any (ncohorts > 0)){
            #----- Make a cohort-level area. ----------------------------------------------#
            areaconow    = rep(areapa,times=ncohorts)
            #------------------------------------------------------------------------------#


            #----- Define the DBH classes. ------------------------------------------------#
            dbhconow     = mymont$DBH
            dbhcut       = cut(dbhconow,breaks=census.obs$dbh.breaks)
            dbhlevs      = levels(dbhcut)
            dbhfac       = match(dbhcut,dbhlevs)
            n.dbh        = length(dbhlevs)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Load the other cohort-level variables of interest.                       #
            #------------------------------------------------------------------------------#
            pftconow       = mymont$PFT
            nplantconow    = mymont$NPLANT
            mortconow      = rowSums(mymont$MMEAN.MORT.RATE)
            ddmortconow    = mymont$MMEAN.MORT.RATE[,2]
            dimortconow    = mortconow - ddmortconow
            recruitconow   = mymont$RECRUIT.DBH
            censtatusconow = mymont$CENSUS.STATUS
            growthconow    = mymont$DLNDBH.DT
            #------------------------------------------------------------------------------#
         }else{
            areaconow       = NA
            dbhconow        = NA
            pftconow        = NA
            nplantconow     = NA
            mortconow       = NA
            ddmortconow     = NA
            dimortconow     = NA
            recruitconow    = NA
            censtatusconow  = NA
            growthconow     = NA
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     The following variable is used to scale "intensive" properties              #
         # (whatever/plant) to "extensive" (whatever/m2).  Sometimes it may be used to     #
         # build weighted averages.                                                        #
         #---------------------------------------------------------------------------------#
         w.nplant = nplantconow * areaconow
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the growth, mortality, and recruitment rates for each PFT, and the     #
         # global rates.  We only use the cohorts that were flagged as 1 or 2 (which means #
         # that their DBH is greater than 10 cm).                                          #
         #---------------------------------------------------------------------------------#
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
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the PFT-level mortality rates.                                      #
            #------------------------------------------------------------------------------#
            if (any(sel.cs2)){

               #---------------------------------------------------------------------------#
               #    Find the weight of each PFT.                                           #
               #---------------------------------------------------------------------------#
               nplantpft     [p] = sum( w.nplant[sel.cs2] )
               #---------------------------------------------------------------------------#

               #---- This is the number of survivors. -------------------------------------#
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
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the PFT-level recruitment rates.                                    #
            #------------------------------------------------------------------------------#
            if (any(sel.dbh) & any(sel.cs2)){
               #---- This is the number of survivors. -------------------------------------#
               population        = sum(w.nplant[sel.dbh])
               established       = sum(w.nplant[sel.cs2])
               recrpft       [p] = log( population / established) / 12.0
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Growth rates are found only for established cohorts.                     #
            #------------------------------------------------------------------------------#
            if (any(sel.cs2)){
               growthpft     [p] = sum( w.nplant[sel.cs2] * growthconow [sel.cs2] )
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Find the mean rates amongst all PFTs.  Because Condit et al. (2006)        #
         # assumed log-normal distribution, we do the same here.  Notice that we don't do  #
         # any weighted mean, and this is because we are looking at the species            #
         # distribution of mortality rates.                                                #
         #---------------------------------------------------------------------------------#
         sel.pft           = nplantpft > 0
         if (any(sel.pft)){
            ts.growth.plot[m] = sum (x=growthpft[sel.pft]) / sum(nplantpft[sel.pft])
            ts.mort.plot  [m] = mean(x=mortpft  [sel.pft])
            ts.dimort.plot[m] = mean(x=dimortpft[sel.pft])
            ts.ddmort.plot[m] = mean(x=ddmortpft[sel.pft])
            ts.recr       [m] = mean(x=recrpft  [sel.pft])
         }else{
            ts.growth.plot[m] = NA
            ts.mort.plot  [m] = NA
            ts.dimort.plot[m] = NA
            ts.ddmort.plot[m] = NA
            ts.recr       [m] = NA
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Build the size (DBH) structure arrays.                                      #
         #---------------------------------------------------------------------------------#
         nplantpftdbh = matrix(0.,nrow=n.dbh,ncol=npft)
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
               sel      = selpft & seldbh & censtatusconow == 2
               if (any(sel)){
                  nplantpftdbh [d,p] = sum( w.nplant [sel])

                  growthpftdbh [d,p] = sum( w.nplant [sel] * growthconow [sel] )

                  #---- This is the number of survivors and living before. ----------------#
                  survivor           = sum( w.nplant[sel] )
                  previous           = sum( w.nplant[sel] * exp(mortconow[sel]) )
                  mortpftdbh   [d,p] = log( previous / survivor )

                  survivor           = sum( w.nplant[sel] )
                  previous           = sum( w.nplant[sel] * exp(dimortconow[sel]) )
                  dimortpftdbh [d,p] = log( previous / survivor )

                  survivor           = sum( w.nplant[sel] )
                  previous           = sum( w.nplant[sel] * exp(ddmortconow[sel]) )
                  ddmortpftdbh [d,p] = log( previous / survivor )
               }#end if
               #---------------------------------------------------------------------------#
            }#end for PFT
            #------------------------------------------------------------------------------#
         }#end for DBH
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Find the mean rates amongst all PFTs.  Because Condit et al. (2006)        #
         # assumed log-normal distribution, we do the same here.  Notice that we don't     #
         # do any weighted mean, and this is because we are looking at the species         #
         # distribution of mortality rates.                                                #
         #---------------------------------------------------------------------------------#
         for (d in 1:n.dbh){
            sel.pft           = nplantpftdbh[d,] > 0
            if (any(sel.pft)){
               ts.growth.size[d,m] = mean(x=growthpftdbh[d,sel.pft])
               ts.mort.size  [d,m] = mean(x=mortpftdbh  [d,sel.pft])
               ts.dimort.size[d,m] = mean(x=dimortpftdbh[d,sel.pft])
               ts.ddmort.size[d,m] = mean(x=ddmortpftdbh[d,sel.pft])
            }else{
               ts.growth.size[d,m] = NA
               ts.mort.size  [d,m] = NA
               ts.dimort.size[d,m] = NA
               ts.ddmort.size[d,m] = NA
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for (m in 1:n.months)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the average rates for the census period using log-normal.                 #
      #------------------------------------------------------------------------------------#
      print(paste("   - Finding the average rates...",sep=""))
      model$cwhen     = census.obs$when
      #----- Plot-level. ------------------------------------------------------------------#
      model$recr        = c(NA,tapply( X     = ts.recr
                                     , INDEX = census.idx
                                     , FUN   = mean
                                     , na.rm = TRUE
                                     ) )
      model$mort.plot   = c(NA,tapply( X     = ts.mort.plot
                                     , INDEX = census.idx
                                     , FUN   = mean
                                     , na.rm = TRUE
                                     ) )
      model$dimort.plot = c(NA,tapply( X     = ts.dimort.plot
                                     , INDEX = census.idx
                                     , FUN   = mean
                                     , na.rm = TRUE
                                     ) )
      model$ddmort.plot = c(NA,tapply( X     = ts.ddmort.plot
                                     , INDEX = census.idx
                                     , FUN   = mean
                                     , na.rm = TRUE
                                     ) )
      model$growth.plot = c(NA,tapply( X     = ts.growth.plot
                                     , INDEX = census.idx
                                     , FUN   = mean
                                     , na.rm = TRUE
                                     ) )
      #----- Size- and plot-level. --------------------------------------------------------#
      empty           = rep(NA,times=n.dbh)
      model$mort.size   = cbind(empty,qapply( X     = ts.mort.size
                                            , DIM   = 2
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$dimort.size = cbind(empty,qapply( X     = ts.dimort.size
                                            , DIM   = 2
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$ddmort.size = cbind(empty,qapply( X     = ts.ddmort.size
                                            , DIM   = 2
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      model$growth.size = cbind(empty,qapply( X     = ts.growth.size
                                            , DIM   = 2
                                            , INDEX = census.idx
                                            , FUN   = mean
                                            , na.rm = TRUE
                                            ) )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     In case vegetation dynamics is turned off, all variables above will be NA.     #
      #------------------------------------------------------------------------------------#
      model$recr        [!is.finite(model$recr       )] = 0.
      model$mort.plot   [!is.finite(model$mort.plot  )] = 0.
      model$dimort.plot [!is.finite(model$dimort.plot)] = 0.
      model$ddmort.plot [!is.finite(model$ddmort.plot)] = 0.
      model$growth.plot [!is.finite(model$growth.plot)] = 0.
      #----- Size- and plot-level. --------------------------------------------------------#
      model$mort.size   [!is.finite(model$mort.size  )] = 0.
      model$dimort.size [!is.finite(model$dimort.size)] = 0.
      model$ddmort.size [!is.finite(model$ddmort.size)] = 0.
      model$growth.size [!is.finite(model$growth.size)] = 0.
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
         this.obs        = 100. * census.obs[[this.vnam.obs]][4,2:n.census]
         this.025        = 100. * census.obs[[this.vnam.obs]][5,2:n.census]
         this.975        = 100. * census.obs[[this.vnam.obs]][6,2:n.census]
         this.mod        = 100. * model     [[this.vnam.ed ]][  2:n.census]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the x axis.                                                           #
         #---------------------------------------------------------------------------------#
         whenplot = pretty.time(when,n=6)
         #---------------------------------------------------------------------------------#

         #----- Find the plot range. ------------------------------------------------------#
         ylim.test = c(this.025,this.975,this.mod)
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
         err.x = c(when,rev(when))
         err.y = c(this.025,rev(this.975))
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
                         ,pointsize=ptsz,paper=paper)
            }#end if
            #------------------------------------------------------------------------------#


            #----- Make plot annotation. --------------------------------------------------#
            letitre = paste(lieu,"\n","Average ",this.desc,sep="")
            lex     = paste("Census time")
            ley     = paste(this.desc,this.unit,sep=" ")
            #------------------------------------------------------------------------------#


            #----- Start the plot. --------------------------------------------------------#
            plot(x=when,y=this.obs,type="n",main=letitre,xlab=lex,ylab=ley,ylim=ylimit
                ,log=plog,xaxt="n",cex.main=cex.main)
            #----- Special, time-friendly X-Axis and grid. --------------------------------#
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){
               abline(v=whenplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
            }#end if
            #----- Confidence interval. ---------------------------------------------------#
            epolygon(x=err.x,y=err.y,col=this.col.obser[1],angle=-45,density=40
                    ,lty="solid",lwd=1.0)
            #----- Observed data. ---------------------------------------------------------#
            points(x=when,y=this.obs,col=this.col.obser[2],lwd=3.0,type="o",pch=16,cex=1.0)
            #----- Modelled data. ---------------------------------------------------------#
            points(x=when,y=this.mod,col=this.col.model,lwd=3.0,type="o",pch=16,cex=1.0)
            #----- Legend. ----------------------------------------------------------------#
            legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                  ,fill   = c(this.col.obser[1],"white")
                  ,border = c(this.col.obser[1],"white")
                  ,angle=-45,density=40,lwd=3.0,col=c(this.col.obser[2],this.col.model)
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
         this.obs        = 100. * census.obs[[this.vnam.obs]][4,,2:n.census]
         this.025        = 100. * census.obs[[this.vnam.obs]][5,,2:n.census]
         this.975        = 100. * census.obs[[this.vnam.obs]][6,,2:n.census]
         this.mod        = 100. * model     [[this.vnam.ed ]][  ,2:n.census]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Find the x axis.                                                           #
         #---------------------------------------------------------------------------------#
         whenplot = pretty.time(when,n=6)
         #---------------------------------------------------------------------------------#

         #----- Find the plot range. ------------------------------------------------------#
         ylim.test = c(this.025,this.975,this.mod)
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
            err.x = c(when,rev(when))
            err.y = c(this.025[d,],rev(this.975[d,]))
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
                            ,pointsize=ptsz,paper=paper)
               }#end if
               #---------------------------------------------------------------------------#


               #----- Make plot annotation. -----------------------------------------------#
               letitre = paste(lieu,"\n","Average ",this.desc
                                        ," - DBH Class: ",dbh.names[d],sep="")
               lex     = paste("Census time")
               ley     = paste(this.desc,this.unit,sep=" ")
               #---------------------------------------------------------------------------#


               #----- Start the plot. -----------------------------------------------------#
               plot(x=when,y=this.obs[d,],type="n",main=letitre,xlab=lex,ylab=ley
                   ,ylim=ylimit,log=plog,xaxt="n",cex.main=cex.main)
               #----- Special, time-friendly X-Axis and grid. -----------------------------#
               axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
               if (plotgrid){
                  abline(v=whenplot$levels,h=axTicks(side=2),col="gray52",lty="solid")
               }#end if
               #----- Confidence interval. ------------------------------------------------#
               epolygon(x=err.x,y=err.y,col=this.col.obser[1],angle=-45,density=40
                       ,lty="solid",lwd=1.0)
               #----- Observed data. ------------------------------------------------------#
               points(x=when,y=this.obs[d,],col=this.col.obser[2],lwd=3.0
                     ,type="o",pch=16,cex=1.0)
               #----- Modelled data. ------------------------------------------------------#
               points(x=when,y=this.mod[d,],col=this.col.model,lwd=3.0
                     ,type="o",pch=16,cex=1.0)
               #----- Legend. -------------------------------------------------------------#
               legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                     ,fill   = c(this.col.obser[1],"white")
                     ,border = c(this.col.obser[1],"white")
                     ,angle=-45,density=40,lwd=3.0,col=c(this.col.obser[2],this.col.model)
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
         this.obs        = 100. * census.obs[[this.vnam.obs]][4,,]
         this.025        = 100. * census.obs[[this.vnam.obs]][5,,]
         this.975        = 100. * census.obs[[this.vnam.obs]][6,,]
         this.mod        = 100. * model     [[this.vnam.ed ]][  ,]
         #---------------------------------------------------------------------------------#


         #----- Find the plot range. ------------------------------------------------------#
         ylim.test = c(this.025,this.975,this.mod)
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
            err.x = c(x.dbh,rev(x.dbh))
            err.y = c(this.025[,y],rev(this.975[,y]))
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
                            ,pointsize=ptsz,paper=paper)
               }#end if
               #---------------------------------------------------------------------------#


               #----- Make plot annotation. -----------------------------------------------#
               letitre = paste(lieu,"\n","Average ",this.desc," - Period: ",period,sep="")
               lex     = paste("Minimum DBH for this class [cm]")
               ley     = paste(this.desc,this.unit,sep=" ")
               #---------------------------------------------------------------------------#


               #----- Start the plot. -----------------------------------------------------#
               plot(x=x.dbh,y=this.obs[,y],type="n",main=letitre,xlab=lex,ylab=ley
                   ,ylim=ylimit,log=plog,cex.main=cex.main)
               #----- Plot the grid. ------------------------------------------------------#
               if (plotgrid) grid(col="gray52",lty="solid")
               #----- Confidence interval. ------------------------------------------------#
               epolygon(x=err.x,y=err.y,col=this.col.obser[1],angle=-45,density=40
                       ,lty="solid",lwd=1.0)
               #----- Observed data. ------------------------------------------------------#
               points(x=x.dbh,y=this.obs[,y],col=this.col.obser[2],lwd=3.0,type="o"
                     ,pch=16,cex=1.0)
               #----- Modelled data. ------------------------------------------------------#
               points(x=x.dbh,y=this.mod[,y],col=this.col.model,lwd=3.0,type="o"
                     ,pch=16,cex=1.0)
               #----- Legend. -------------------------------------------------------------#
               legend(x="topleft",inset=0.01,legend=c("Observation","Model")
                     ,fill   = c(this.col.obser[1],"white")
                     ,border = c(this.col.obser[1],"white")
                     ,angle=-45,density=40,lwd=3.0,col=c(this.col.obser[2],this.col.model)
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
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#

#q("no")
