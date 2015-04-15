#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
here    = getwd()                                #   Current directory
srcdir  = "/n/home00/mlongo/util/Rsc"            #   Script directory
outroot = file.path(here,"radprof_comp")         #   Output directory

#------ Plot options. ---------------------------------------------------------------------#
outform        = c("pdf")              # Formats for output file.  Supported formats are:
                                       #   - "X11" - for printing on screen
                                       #   - "eps" - for postscript printing
                                       #   - "png" - for PNG printing
                                       #   - "pdf" - for PDF printing
depth          = 96                    # PNG resolution, in pixels per inch
paper          = "letter"              # Paper size, to define the plot shape
wpaper         = "legal"               # Wide paper size, to define the plot shape
ptsz           = 18                    # Font size.
ibackground    = 0                     # Make figures compatible to which background?
                                       # 0 -- white
                                       # 1 -- black
                                       # 2 -- dark grey
#------------------------------------------------------------------------------------------#




#----- Scale shortwave variables? ---------------------------------------------------------#
scale.sw = FALSE
#------------------------------------------------------------------------------------------#


#----- Profile options. -------------------------------------------------------------------#
my.hours   = sequence(24)-1
prof.hours = c(0,3,6,9,12,15,18,21)                  # Hours to try
prof.zen   = c(0,5,20,30,40,50,60,70,75,80,84,88,92) # Zenith angle breaks
#------------------------------------------------------------------------------------------#


#----- Default break for LAI classes. -----------------------------------------------------#
lai.defbrk  = c(  1.5,  2.0,  2.5,  3.0,  3.5,  4.5,  5.5,  7.0)
add.others  = FALSE
report      = list()
#---------------------------------------   1.5   2.0   2.5   3.0   3.5   4.5   5.5   7.0  -#
report[[1]] = list( solar   = rbind( min = c(  NA ,  NA , 22.0,  8.9,  5.4,  NA ,  NA )
                                   , max = c(  NA ,  NA , 45.4, 11.3,  6.7,  NA ,  NA )
                                   )#end rbind
                  , par     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , nir     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , colour  = "steelblue1"
                  , angle   = 20
                  , density = 40
                  , lwd     = 2
                  , desc    = "Torquebiau and Walter (1987)"
                  )#end report
#---------------------------------------   1.5   2.0   2.5   3.0   3.5   4.5   5.5   7.0  -#
report[[2]] = list( solar   = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , par     = rbind( min = c(  NA ,  NA , 10.8,  2.6, 0.35,  NA ,  NA )
                                   , max = c(  NA ,  NA , 22.2,  3.3, 1.04,  NA ,  NA )
                                   )#end rbind
                  , nir     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , colour  = "mediumpurple1"
                  , angle   = 70
                  , density = 20
                  , lwd     = 2
                  , desc    = "Torquebiau (1988)"
                  )#end report
#---------------------------------------   1.5   2.0   2.5   3.0   3.5   4.5   5.5   7.0  -#
report[[3]] = list( solar   = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA , 0.89)
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA , 1.61)
                                   )#end rbind
                  , par     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , nir     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , colour  = "sandybrown"
                  , angle   = -20
                  , density =  40
                  , lwd     =   2
                  , desc    = "Shuttleworth (1984)"
                  )#end report
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     List with all sites.  List elements are:                                             #
#     - iata: site code                                                                    #
#     - desc: full name                                                                    #
#     - pch : symbol for this site (currently ignored)                                     #
#     - col : default colour for this site (currently ignored)                             #
#     - year: year to process                                                              #
#                                                                                          #
# The use.sites key allows to subset.  It can be either the number or the iata code.  For  #
#     convenience, you can also set it to TRUE, which will load all sites (FALSE will skip #
#     all of them, kind of useless...)                                                     #
#------------------------------------------------------------------------------------------#
sites      = list()
sites[[1]] = list(iata="gyf",desc="Paracou"       ,lai=6.0,pch= 2,col="#520485",year=2004)
sites[[2]] = list(iata="s67",desc="Santarem km 67",lai=4.0,pch= 5,col="#46FF32",year=2009)
sites[[3]] = list(iata="s83",desc="Santarem km 83",lai=5.0,pch= 9,col="#FF5700",year=2001)
sites[[4]] = list(iata="pdg",desc="Pe-de-Gigante" ,lai=3.5,pch=13,col="#A00014",year=2003)
sites[[5]] = list(iata="rja",desc="Rebio Jaru"    ,lai=4.0,pch= 1,col="#006715",year=2001)
sites[[6]] = list(iata="m34",desc="Manaus K34"    ,lai=5.8,pch= 6,col="#0742C3",year=2004)
sites[[7]] = list(iata="pnz",desc="Petrolina"     ,lai=2.0,pch= 4,col="#B49ED2",year=2004)
sites[[8]] = list(iata="ban",desc="Bananal"       ,lai=3.0,pch= 8,col="#F5C858",year=2005)
sites[[9]] = list(iata="dcm",desc="Santarem km 67",lai=4.0,pch= 5,col="#46FF32",year=2009)
use.sites  = "dcm" # "rja"
imetrad    = paste0("imetrad",sprintf("%2.2i",c(2,3,4)))[3]
#------------------------------------------------------------------------------------------#





#----- Info on hourly data. ---------------------------------------------------------------#
reload.profile  = TRUE
rdata.path      = file.path(here,"RData_profile")
rdata.suffix    = paste0("radprof_",imetrad,"_ed22.RData")
finished.suffix = paste0("radprof_",imetrad,"_ed22.txt")
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.summary  = TRUE
rdata.summary   = paste0("summary_radprof_",imetrad,"_ed22.RData")
#------------------------------------------------------------------------------------------#



#----- Flags to control whether to plot site-specific and multi-site comparison. ----------#
plot.site    = c(FALSE,TRUE)[2]
plot.patch   = c(FALSE,TRUE)[1]
plot.month   = c(FALSE,TRUE)[2]
plot.refonly = c(FALSE,TRUE)[2]
dlai        = 0.25
ci.estimate = "six.summary"    # How to define the statistics?
                               # "boot.summary" -- Use bootstrap to get CI
                               # "six.summary"  -- Use t distribution to get CI
n.boot      = 1000
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



#----- Load profiles from ABRACOS. --------------------------------------------------------#
load(file.path(srcdir,"ABRACOS.light.RData"))
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     List all simulations.                                                                #
#------------------------------------------------------------------------------------------#
simul       = list()
simul[[ 1]] = list( suff           = paste("labpar600",imetrad,sep="_")
                  , desc           = "LA600"
                  , colour         = "#0059B3"
                  , bgcol          = "#537CA6"
                  , mult.down      = 1.
                  , mult.beam.down = 1.
                  , mult.diff.down = 0.97
                  , mult.up        = 1.
                  , angle          = -30
                  , density        = 25
                  )#end list
simul[[ 2]] = list( suff           = paste("labpar850",imetrad,sep="_")
                  , desc           = "LA850"
                  , colour         = "#4D0099"
                  , bgcol          = "#8659B3"
                  , mult.down      = 1.
                  , mult.beam.down = 1.
                  , mult.diff.down = 1.03
                  , mult.up        = 0.985
                  , angle          = 30
                  , density        = 25
                  )#end list
simul[[ 3]] = list( suff           = paste("labpar925",imetrad,sep="_")
                  , desc           = "LA925"
                  , colour         = "#B3B300"
                  , bgcol          = "#A6A653"
                  , mult.down      = 1.03
                  , mult.beam.down = 1.03
                  , mult.diff.down = 1.
                  , mult.up        = 1.
                  , angle          = -60
                  , density        = 25
                  )#end list
simul[[ 4]] = list( suff           = paste("labpar980",imetrad,sep="_")
                  , desc           = "LA980"
                  , colour         = "#446600"
                  , bgcol          = "#AACC66"
                  , mult.down      = 0.97
                  , mult.beam.down = 0.97
                  , mult.diff.down = 1.
                  , mult.up        = 1.015
                  , angle          = 60
                  , density        = 25
                  )#end list
#------------------------------------------------------------------------------------------#


#----- List of bands. ---------------------------------------------------------------------#
band = list()
band[[1]] =  list( key  = "sol"
                 , ed2  = "RSHORT"
                 , desc = "Total Solar"
                 , unit = "wom2"
                 , col  = "#0742C3"
                 )#end list
band[[2]] =  list( key  = "par"
                 , ed2  = "PAR"
                 , desc = "Photosynthetically Active"
                 , unit = "umolom2os"
                 , col  = "#46FF32"
                 )#end list
band[[3]] =  list( key  = "nir"
                 , ed2  = "NIR"
                 , desc = "Near Infrared"
                 , unit = "wom2"
                 , col  = "#FF5700"
                 )#end list
band[[4]] =  list( key  = "tir"
                 , ed2  = "RLONG"
                 , desc = "Thermal Infrared"
                 , unit = "wom2"
                 , col  = "#70000E"
                 )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Convert lists to a data frames.  For sites only, we may skip several of them using  #
# use.sites.                                                                               #
#------------------------------------------------------------------------------------------#
sites = list.2.data.frame(sites)
if (is.logical(use.sites)){
   use.sites = rep(use.sites,times=nrow(sites))
}else if (is.character(use.sites)){
   use.sites = match(tolower(use.sites),sites$iata)
   use.sites = use.sites[is.finite(use.sites)]
}#end if
sites = sites[use.sites,]
simul = list.2.data.frame(simul)
band  = list.2.data.frame(band )
#------------------------------------------------------------------------------------------#



#----- Number of sites and simulations. ---------------------------------------------------#
nsites  = nrow(sites)
nsimul  = nrow(simul)
nband   = nrow(band )
nprof   = length(prof.hours)
nhours  = length(my.hours)
nreport = length(report)
#------------------------------------------------------------------------------------------#



#----- Number of output formats. ----------------------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size  = plotsize(proje=FALSE,paper=paper )
wsize = plotsize(proje=FALSE,paper=wpaper)
tsize = plotsize(proje=FALSE,stdheight=8.5*6/5,stdwidth=11)
#------------------------------------------------------------------------------------------#




#----- Define sub-plot instructions for simulation comparison. ----------------------------#
lo.simul = pretty.box(n=nsimul)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
metroot = file.path(outroot,imetrad)
if (! file.exists(rdata.path)) dir.create(rdata.path)
if (! file.exists(outroot   )) dir.create(outroot   )
if (! file.exists(metroot   )) dir.create(metroot   )
#------------------------------------------------------------------------------------------#


#----- Summary file. ----------------------------------------------------------------------#
summary.fullname = file.path(rdata.path,rdata.summary)
#------------------------------------------------------------------------------------------#



#----- Zenith angle bins. -----------------------------------------------------------------#
nzb       = length(prof.zen)
zen.key   = paste("zen_",sprintf("%2.2i",prof.zen[-nzb]),"-"
                        ,sprintf("%2.2i",prof.zen[-1]),sep="")
zen.desc  = paste("Mean zenith angle ",sprintf("%.1f",mid.points(prof.zen)),sep="")
nzen.prof = length(zen.key)
#------------------------------------------------------------------------------------------#


#----- These files are quite heavy, so we load them only if they are necessary. -----------#
if ( (! file.exists(summary.fullname)) || plot.site){
   #---------------------------------------------------------------------------------------#
   #      Create all output directories, separated by format.                              #
   #---------------------------------------------------------------------------------------#
   eft        = list()
   loop.sites = integer(length=0)
   for (p in sequence(nsites)){
      #----- Retrieve site sinformation. --------------------------------------------------#
      iata     = sites$iata[p]
      ip       = match(iata,poilist$iata)
      yeara    = sites$year[p]
      yearz    = sites$year[p]
      lon      = poilist$lon     [ip]
      lat      = poilist$lat     [ip]
      desc     = poilist$longname[ip]
      short    = poilist$short   [ip]
      #------------------------------------------------------------------------------------#



      #----- Find file name. --------------------------------------------------------------#
      rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
      if (file.exists(rdata.iata)){
         #----- Reload data and copy to the general list. ---------------------------------#
         cat(" + Loading data from ",basename(rdata.iata),"...","\n",sep="")
         dummy       = load(file=rdata.iata)
         eft[[iata]] = model
         rm(model)
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #      Find all days that we will grab data.                                      #
         #---------------------------------------------------------------------------------#
         cat(" - Generating profiles for ",desc,"...","\n")
         alldays   = chron(seq( from = as.numeric(chron(paste( 1, 1,yeara,sep="/")) )
                              , to   = as.numeric(chron(paste(12,31,yearz,sep="/")) )
                              , by   = 1
                              )#end seq
                          )#end chron
         ndates    = length(alldays)
         when      = chron( dates = rep(paste(alldays),each=nhours)
                          , times = paste(rep(my.hours,times=ndates),0,0,sep=":")
                          )#end chron
         w.year    = numyears (when)
         w.month   = nummonths(when)
         w.day     = numdays  (when)
         w.hour    = hours    (when)
         w.tomonth = chron(paste(w.month,1,w.year,sep="/"))
         ntimes    = sort(unique(w.tomonth))


         zenith   = mapply( FUN      = function(bef,now,...){
                                          zz  = ed.zen(when=c(bef,now),...)
                                          ans = c( zen   = zz$zen  [2]
                                                 , day   = zz$day  [2]
                                                 , night = zz$night[2]
                                                 )#end c
                                          return(ans)
                                       }#end function
                          , bef      = when - 1/24
                          , now      = when
                          , MoreArgs = list( lon       = lon
                                           , lat       = lat
                                           , ed21      = TRUE
                                           , zeronight = FALSE
                                           , meanval   = TRUE
                                           , imetavg   = 1
                                           , nmean     = 15
                                           )#end list
                          , SIMPLIFY = TRUE
                          )#end mapply

         zen   = tapply( X     = zenith["zen",]
                       , INDEX = list(w.tomonth,w.hour)
                       , FUN   = mean
                       , na.rm = TRUE
                       )#end tapply
         day   = tapply( X     = zenith["day",]
                       , INDEX = list(w.tomonth,w.hour)
                       , FUN   = commonest
                       , na.rm = TRUE
                       )#end tapply
         night = tapply( X     = zenith["night",]
                       , INDEX = list(w.tomonth,w.hour)
                       , FUN   = commonest
                       , na.rm = TRUE
                       )#end tapply

         model    = list( iata         = iata
                        , short        = short
                        , longname     = desc
                        , lon          = lon
                        , lat          = lat
                        , yeara        = yeara
                        , yearz        = yearz
                        , tresume      = 1
                        , tomonth      = sort(unique(w.tomonth))
                        , month        = nummonths(sort(unique(w.tomonth)))
                        , year         = numyears (sort(unique(w.tomonth)))
                        , hour         = sort(unique(w.hour))
                        , ntimes       = length(unique(w.tomonth))
                        , ndcycle      = length(unique(w.hour))
                        , zen          = zen
                        , diel         = round( 1 + day + night) + 0.*night
                        , gloom        = matrix( data = FALSE
                                               , ncol = ncol(night)
                                               , nrow = nrow(night)
                                               )#end matrix
                        )#end list
         #---------------------------------------------------------------------------------#

         eft[[iata]] = model

      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#





   #=======================================================================================#
   #=======================================================================================#
   #      Loop over sites.                                                                 #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(nsites)){

      #------ Grab data for this site. ----------------------------------------------------#
      iata      = sites$iata[p]
      model     = eft[[iata]]
      yeara     = model$yeara
      yearz     = model$yearz
      lon       = model$lon
      lat       = model$lat
      longname  = model$longname
      short     = model$short
      tresume   = model$tresume
      ntimes    = model$ntimes
      ndcycle   = model$ndcycle
      tomonth   = model$tomonth
      #------------------------------------------------------------------------------------#


      #----- File name for this site. -----------------------------------------------------#
      rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
      #------------------------------------------------------------------------------------#

      #---- Get the times that must be read. ----------------------------------------------#
      if (tresume > ntimes){
         loop.times = numeric(0)
      }else{
         loop.times = seq(from=tresume,to=ntimes,by=1)
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over dates.                                                               #
      #------------------------------------------------------------------------------------#
      for (w in loop.times){


         #----- Grab date/time information on the data set. -------------------------------#
         now   = tomonth[w]
         yyyy = sprintf("%4.4i",numyears (now))
         mm   = sprintf("%2.2i",nummonths(now))
         dd   = sprintf("%2.2i",numdays  (now))
         hh   = sprintf("%2.2i",hours    (now))
         cat("   * ",paste(now),"\n",sep="")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#



         #----- Build file suffix. --------------------------------------------------------#
         tsuff = paste(yyyy,"-",mm,"-00-000000-g01.h5",sep="")
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Loop over simulations.                                                       #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nsimul)){

            #----- Build the file name. ---------------------------------------------------#
            simplace = paste("t",iata,"_",simul$suff[s],sep="")
            filenow  = paste(simplace,"Q",tsuff,sep="-")
            #------------------------------------------------------------------------------#


            #----- Make file name as well as the compressed names. ------------------------#
            fullfile     = file.path(here,simplace,"analy",filenow)
            fullfile.gz  = paste(fullfile,"gz" ,sep=".")
            fullfile.bz2 = paste(fullfile,"bz2",sep=".")
            #------------------------------------------------------------------------------#



            #----- Open the file. ---------------------------------------------------------#
            if (file.exists(fullfile)){
               mymont    = hdf5load(file=fullfile,load=FALSE,tidy=TRUE,verbosity=0)
            }else if (file.exists(fullfile.gz )){
               temp.file = file.path(tempdir(),basename(fullfile))
               dummy     = gunzip(filename=fullfile.gz,destname=temp.file,remove=FALSE)
               mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else if (file.exists(fullfile.bz2)){
               temp.file = file.path(tempdir(),basename(fullfile))
               dummy     = bunzip2(filename=fullfile.bz2,destname=temp.file,remove=FALSE)
               mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else{
               stop ("File not found")
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Append some auxiliary variables to mymont.                              #
            #------------------------------------------------------------------------------#
            mymont$SIPA.N                   = length(mymont$PACO.ID)
            mymont$QMEAN.ATM.NIR.PY         = ( mymont$QMEAN.ATM.RSHORT.PY
                                              - mymont$QMEAN.ATM.PAR.PY         )
            mymont$QMEAN.ATM.NIR.DIFF.PY    = ( mymont$QMEAN.ATM.RSHORT.DIFF.PY
                                              - mymont$QMEAN.ATM.PAR.DIFF.PY    )
            mymont$QMEAN.ATM.RSHORT.BEAM.PY = ( mymont$QMEAN.ATM.RSHORT.PY
                                              - mymont$QMEAN.ATM.RSHORT.DIFF.PY )
            mymont$QMEAN.ATM.PAR.BEAM.PY    = ( mymont$QMEAN.ATM.PAR.PY
                                              - mymont$QMEAN.ATM.PAR.DIFF.PY    )
            mymont$QMEAN.ATM.NIR.BEAM.PY    = ( mymont$QMEAN.ATM.NIR.PY
                                              - mymont$QMEAN.ATM.NIR.DIFF.PY    )
            mymont$QMEAN.ATM.RLONG.DIFF.PY  = mymont$QMEAN.ATM.RLONG.PY
            mymont$QMEAN.ATM.RLONG.BEAM.PY  = 0. * mymont$QMEAN.ATM.RLONG.PY
            mymont$QMEAN.ATM.NIR.SI         = ( mymont$QMEAN.ATM.RSHORT.SI
                                              - mymont$QMEAN.ATM.PAR.SI         )
            mymont$QMEAN.ATM.NIR.DIFF.SI    = ( mymont$QMEAN.ATM.RSHORT.DIFF.SI
                                              - mymont$QMEAN.ATM.PAR.DIFF.SI    )
            mymont$QMEAN.ATM.RSHORT.BEAM.SI = ( mymont$QMEAN.ATM.RSHORT.SI
                                              - mymont$QMEAN.ATM.RSHORT.DIFF.SI )
            mymont$QMEAN.ATM.PAR.BEAM.SI    = ( mymont$QMEAN.ATM.PAR.SI
                                              - mymont$QMEAN.ATM.PAR.DIFF.SI    )
            mymont$QMEAN.ATM.NIR.BEAM.SI    = ( mymont$QMEAN.ATM.NIR.SI
                                              - mymont$QMEAN.ATM.NIR.DIFF.SI    )
            mymont$QMEAN.ATM.RLONG.DIFF.SI  = mymont$QMEAN.ATM.RLONG.SI
            mymont$QMEAN.ATM.RLONG.BEAM.SI  = 0. * mymont$QMEAN.ATM.RLONG.SI
            mymont$QMEAN.RSHORTUP.PA        = ( mymont$QMEAN.PARUP.PA 
                                              + mymont$QMEAN.NIRUP.PA )
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the radiation at the top of the atmosphere a patch variables.       #
            #------------------------------------------------------------------------------#
            mymont$QMEAN.PAR.BEAM.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.PAR.BEAM.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.PAR.DIFF.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.PAR.DIFF.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.NIR.BEAM.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.NIR.BEAM.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.NIR.DIFF.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.NIR.DIFF.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RSHORT.BEAM.PA = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RSHORT.BEAM.SI
                                                        , each = SIPA.N
                                                        )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RSHORT.DIFF.PA = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RSHORT.DIFF.SI
                                                        , each = SIPA.N
                                                        )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RLONG.BEAM.PA  = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RLONG.BEAM.SI
                                                         , each = SIPA.N
                                                         )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RLONG.DIFF.PA  = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RLONG.DIFF.SI
                                                         , each = SIPA.N
                                                         )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Correct profiles.                                                       #
            #------------------------------------------------------------------------------#
            mymont$QMEAN.PAR.BEAM.CO    =   mymont$QMEAN.RAD.PROFILE[,, 1]
            mymont$QMEAN.PAR.DIFF.CO    =   mymont$QMEAN.RAD.PROFILE[,, 3]
            mymont$QMEAN.PARUP.CO       = ( mymont$QMEAN.RAD.PROFILE[,, 2]
                                          + mymont$QMEAN.RAD.PROFILE[,, 4] )
            mymont$QMEAN.NIR.BEAM.CO    =   mymont$QMEAN.RAD.PROFILE[,, 5]
            mymont$QMEAN.NIR.DIFF.CO    =   mymont$QMEAN.RAD.PROFILE[,, 7]
            mymont$QMEAN.NIRUP.CO       = ( mymont$QMEAN.RAD.PROFILE[,, 6]
                                          + mymont$QMEAN.RAD.PROFILE[,, 8] )
            mymont$QMEAN.RLONG.BEAM.CO  =   mymont$QMEAN.RAD.PROFILE[,, 9] * 0.
            mymont$QMEAN.RLONG.DIFF.CO  =   mymont$QMEAN.RAD.PROFILE[,, 9]
            mymont$QMEAN.RLONGUP.CO     =   mymont$QMEAN.RAD.PROFILE[,,10]
            mymont$QMEAN.RSHORT.BEAM.CO = ( mymont$QMEAN.PAR.BEAM.CO
                                          + mymont$QMEAN.NIR.BEAM.CO
                                          )#end
            mymont$QMEAN.RSHORT.DIFF.CO = ( mymont$QMEAN.PAR.DIFF.CO
                                          + mymont$QMEAN.NIR.DIFF.CO
                                          )#end
            mymont$QMEAN.RSHORTUP.CO    = ( mymont$QMEAN.PARUP.CO
                                          + mymont$QMEAN.NIRUP.CO
                                          )#end
            mymont$QMEAN.NIR.L.CO       = ( mymont$QMEAN.RSHORT.L.CO
                                          - mymont$QMEAN.PAR.L.CO
                                          )#end
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Load the patch and cohort information.                                  #
            #------------------------------------------------------------------------------#
            if (w == 1 && s == 1){
               nsites    = mymont$PYSI.N
               npatches  = mymont$SIPA.N
               ncohorts  = mymont$PACO.N
               a.ico     = mymont$PACO.ID
               z.ico     = a.ico + ncohorts - 1

               isi       = rep(sequence(nsites),times=npatches)
               ipa       = sequence(npatches)
               ipaco     = rep(ipa,times=ncohorts  )
               ppaco     = rep(ipa,times=ncohorts+1)
               ico       = sequence(ncohorts)
               pco       = sequence(ncohorts+1)-1

               nlayer    = sum(ncohorts+1)

               #----- Initialise data regarding this polygon. -----------------------------#
               empty.ai    = array (data= NA,dim=c(nsimul,ntimes        ,nlayer))
               empty.atm   = array (data= NA,dim=c(nsimul,ntimes,ndcycle       ))
               empty.prof  = array (data= NA,dim=c(nsimul,ntimes,ndcycle,nlayer))
               empty.list  = list  ( down.beam = empty.prof
                                   , down.diff = empty.prof
                                   , down.tot  = empty.prof
                                   , up.tot    = empty.prof
                                   , abs.lyr   = empty.prof
                                   , abs.tot   = empty.prof
                                   , abs.leaf  = empty.prof
                                   , abs.wood  = empty.prof
                                   , atm.tot   = empty.atm
                                   , atm.beam  = empty.atm
                                   , atm.diff  = empty.atm
                                   )#end list
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Append initial conditions to this polygon.                           #
               #---------------------------------------------------------------------------#
               model = modifyList( x   = model
                                 , val = list( npatches = npatches
                                             , ncohorts = ncohorts
                                             , desert   = sum(ncohorts) == 0
                                             , nlayer   = nlayer
                                             , ipa      = ipa
                                             , ipaco    = ipaco
                                             , ppaco    = ppaco
                                             , ico      = ico
                                             , pco      = pco
                                             , lai      = empty.ai
                                             , wai      = empty.ai
                                             , tai      = empty.ai
                                             , par      = empty.list
                                             , nir      = empty.list
                                             , sol      = empty.list
                                             , tir      = empty.list
                                             )#end list
                                 )#end modifyList
               #---------------------------------------------------------------------------#
            }#end if (w == 1 && s == 1)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Retrieve some useful data from this run.                                 #
            #------------------------------------------------------------------------------#
            npatches = model$npatches
            ipa      = model$ipa
            ipaco    = model$ipaco
            ico      = model$ico
            ppaco    = model$ppa
            pco      = model$pco
            desert   = model$desert
            #------------------------------------------------------------------------------#



            #------ Variables to fill in constant boundary conditions. --------------------#
            zero     = rep(0   ,times=npatches)
            veritas  = rep(TRUE,times=npatches)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the cumulative LAI profile.                                         #
            #------------------------------------------------------------------------------#
            if (desert){
               lai.prof = zero
               wai.prof = zero
               tai.prof = zero
               use.lyr  = veritas
            }else{
               lai.prof = append.patch(ipaco=ipaco,xpa=zero,xco=mymont$MMEAN.LAI.CO)
               wai.prof = append.patch(ipaco=ipaco,xpa=zero,xco=mymont$WAI.CO)
               tai.prof = lai.prof + wai.prof
               use.lyr  = append.patch( ipaco = ipaco
                                      , xpa   = veritas
                                      , xco   = mymont$MMEAN.RLONG.L.CO != 0.
                                      )#end append.patch
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Append layer characteristics for this time.                              #
            #------------------------------------------------------------------------------#
            model$lai      [s,w,] = lai.prof
            model$wai      [s,w,] = wai.prof
            model$tai      [s,w,] = tai.prof
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over radiation bands.                                               #
            #------------------------------------------------------------------------------#
            for (b in sequence(nband)){
               bnd  = band$key [b]
               ed2  = band$ed2 [b]



               #----- Get the components. -------------------------------------------------#
               down.beam.top = mymont[[paste("QMEAN.",ed2,".BEAM.PA",sep="")]]
               down.diff.top = mymont[[paste("QMEAN.",ed2,".DIFF.PA",sep="")]]
               up.tot.top    = mymont[[paste("QMEAN.",ed2,"UP.PA"   ,sep="")]]
               down.beam.coh = mymont[[paste("QMEAN.",ed2,".BEAM.CO",sep="")]]
               down.diff.coh = mymont[[paste("QMEAN.",ed2,".DIFF.CO",sep="")]]
               up.tot.coh    = mymont[[paste("QMEAN.",ed2,"UP.CO"   ,sep="")]]
               abs.leaf.coh  = mymont[[paste("QMEAN.",ed2,".L.CO"   ,sep="")]]

               atm.beam      = mymont[[paste("QMEAN.ATM.",ed2,".BEAM.PY",sep="")]]
               atm.diff      = mymont[[paste("QMEAN.ATM.",ed2,".DIFF.PY",sep="")]]
               atm.tot       = atm.beam + atm.diff
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #    Append data by time of day.                                            #
               #---------------------------------------------------------------------------#
               for (h in sequence(ndcycle)){

                  #------------------------------------------------------------------------#
                  #     Add top boundary condition to each cohort.                         #
                  #------------------------------------------------------------------------#
                  if (desert){
                     #----- Desert run, no cohorts. ---------------------------------------#
                     down.beam = down.beam.top[,h]
                     down.diff = down.diff.top[,h]
                     up.tot    = up.top       [,h]
                     abs.leaf  = zero
                     #---------------------------------------------------------------------#
                  }else{
                     #----- Append the patch-level top boundary condition. ----------------#
                     down.beam = append.patch( ipaco = ipaco
                                             , xpa   = down.beam.top[,h]
                                             , xco   = down.beam.coh[,h]
                                             )#end append.patch
                     down.diff = append.patch( ipaco = ipaco
                                             , xpa   = down.diff.top[,h]
                                             , xco   = down.diff.coh[,h]
                                             )#end append.patch
                     up.tot    = append.patch( ipaco = ipaco
                                             , xpa   = up.tot.top   [,h]
                                             , xco   = up.tot.coh   [,h]
                                             )#end append.patch
                     abs.leaf  = append.patch( ipaco = ipaco
                                             , xpa   = zero
                                             , xco   = abs.leaf.coh [,h]
                                             )#end append.patch
                     #---------------------------------------------------------------------#



                     #----- Fill in cohorts that haven't been resolved. -------------------#
                     down.beam = fill.unresolved(ipaco=ppaco,x=down.beam,use=use.lyr)
                     down.diff = fill.unresolved(ipaco=ppaco,x=down.diff,use=use.lyr)
                     up.tot    = fill.unresolved(ipaco=ppaco,x=up.tot   ,use=use.lyr)
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Find the derived quantities. -------------------------------------#
                  down.tot  = down.beam + down.diff
                  abs.lyr   = layer.absorption(ipaco=ppaco,down=down.tot,up=up.tot)
                  abs.wood  = abs.lyr - abs.leaf
                  abs.tot   = unlist(tapply(X=abs.lyr ,INDEX=ppaco,FUN=cumsum))
                  abs.leaf  = unlist(tapply(X=abs.leaf,INDEX=ppaco,FUN=cumsum))
                  abs.wood  = unlist(tapply(X=abs.wood,INDEX=ppaco,FUN=cumsum))
                  #------------------------------------------------------------------------#



                  #----- Update profile. --------------------------------------------------#
                  model[[bnd]]$down.beam [s,w,h,] = down.beam
                  model[[bnd]]$down.diff [s,w,h,] = down.diff
                  model[[bnd]]$down.tot  [s,w,h,] = down.tot
                  model[[bnd]]$up.tot    [s,w,h,] = up.tot
                  model[[bnd]]$abs.lyr   [s,w,h,] = abs.lyr
                  model[[bnd]]$abs.tot   [s,w,h,] = abs.tot
                  model[[bnd]]$abs.leaf  [s,w,h,] = abs.leaf
                  model[[bnd]]$abs.wood  [s,w,h,] = abs.wood
                  #------------------------------------------------------------------------#
               }#end for (h in sequence(ndcycle))
               #---------------------------------------------------------------------------#



               #----- Update ToC radiation. -----------------------------------------------#
               model[[bnd]]$atm.beam  [s,w,] =  atm.beam
               model[[bnd]]$atm.diff  [s,w,] =  atm.diff
               model[[bnd]]$atm.tot   [s,w,] =  atm.tot
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#

         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (w in sequence(nwhen))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#
      cat("     > Prepare to save...","\n",sep="")
      if (tresume <= ntimes){
         #---------------------------------------------------------------------------------#
         #     Now we find which patch has the most, the least and the median LAI.         #
         #---------------------------------------------------------------------------------#
         model$lai.pa   = qapply(X=model$lai,INDEX =model$ppaco,DIM=3,FUN=sum ,na.rm=TRUE)
         lai.bar        = apply (X=model$lai.pa,MARGIN=3,FUN=mean,na.rm=TRUE)
         model$lai.plot = c( max = which.max(lai.bar)
                           , mid = which.closest(x=median(lai.bar),A=lai.bar)
                           , min = which.min(lai.bar)
                           )#end c
         #---------------------------------------------------------------------------------#



         #----- Reload data and copy to the general list. ---------------------------------#
         rdata.iata     = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
         rdata.finished = file.path(rdata.path
                                   ,paste("loaded",iata,finished.suffix,sep="_"))
         cat(" + Saving data to ",basename(rdata.iata),"...","\n",sep="")
         model$tresume = ntimes+1
         dummy         = save(model,file=rdata.iata)
         dummy         = write(x=paste(Sys.time()),file=rdata.finished)
         eft[[iata]]   = model
         rm(model)
         #---------------------------------------------------------------------------------#
      }#end if (tresume <= nwhen)
   }#end for (p in sequence(nsites))
   #=======================================================================================#
   #=======================================================================================#
}#end if ( ( (! file.exists(summary.fullname)) && plot.patch ) || plot.site )
#==========================================================================================#
#==========================================================================================#








#==========================================================================================#
#==========================================================================================#
#     Patch-level summary for all sites.                                                   #
#------------------------------------------------------------------------------------------#
if (! file.exists(summary.fullname) && plot.patch){
   summ = list()


   for (p in sequence(nsites)){
      #----- Retrieve site information. ---------------------------------------------------#
      iata     = sites$iata[p]
      model    = eft[[iata]]
      yeara    = model$yeara
      yearz    = model$yearz
      lon      = model$lon
      lat      = model$lat
      desc     = model$longname
      cat(" - Generating patch-level summary for ",desc,"...","\n")
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Load LAI/WAI/TAI profiles.                                                     #
      #------------------------------------------------------------------------------------#
      lai     = qapply(X=model$lai,INDEX=model$month,DIM=2,FUN=mean)
      wai     = qapply(X=model$wai,INDEX=model$month,DIM=2,FUN=mean)
      tai     = qapply(X=model$tai,INDEX=model$month,DIM=2,FUN=mean)
      tot.lai = qapply(X=model$lai,INDEX=model$ppaco,DIM=3,FUN=sum )
      tot.wai = qapply(X=model$wai,INDEX=model$ppaco,DIM=3,FUN=sum )
      tot.tai = qapply(X=model$tai,INDEX=model$ppaco,DIM=3,FUN=sum )
      #------------------------------------------------------------------------------------#



      #----- Initialise arrays that will hold the mean diel. ------------------------------#
      npatches = model$npatches
      empty    = array( data     = NA
                      , dim      = c(nband,nsimul,12,ndcycle,npatches)
                      , dimnames = list( band$key
                                       , simul$suff
                                       , month.abb
                                       , model$w.hour
                                       , sequence(npatches)
                                       )#end list
                      )#end array
      flux     = list( down.beam = empty
                     , down.diff = empty
                     , down.tot  = empty
                     , up.tot    = empty
                     )#end list
      #------------------------------------------------------------------------------------#




      #----- Band loop. -------------------------------------------------------------------#
      for (b in sequence(nband)){
         bnd      = band$key [b]
         bnd.desc = band$desc[b]
         prof     = model[[bnd]]


         #----- Find the annual averages for each cohort. ---------------------------------#
         down.beam = qapply(X=prof$down.beam,INDEX=model$month,DIM=2,FUN=mean)
         down.diff = qapply(X=prof$down.diff,INDEX=model$month,DIM=2,FUN=mean)
         down.tot  = qapply(X=prof$down.tot ,INDEX=model$month,DIM=2,FUN=mean)
         up.tot    = qapply(X=prof$up.tot   ,INDEX=model$month,DIM=2,FUN=mean)
         #---------------------------------------------------------------------------------#

         #----- Get only the bottom of the cohort. ----------------------------------------#
         flux$down.beam[b,,,] = qapply(X=down.beam,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         flux$down.diff[b,,,] = qapply(X=down.diff,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         flux$down.tot [b,,,] = qapply(X=down.tot ,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         flux$up.tot   [b,,,] = qapply(X=up.tot   ,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         #---------------------------------------------------------------------------------#
      }#end for (b in sequence(nband))
      #------------------------------------------------------------------------------------#



      #----- Save the summary for this site. ----------------------------------------------#
      summ[[iata]] = list( lai = tot.lai, wai = tot.wai, tai = tot.tai, flux = flux )
      #------------------------------------------------------------------------------------#

   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#

   cat(" - Saving patch-level summary...","\n")
   dummy = save(summ,file=summary.fullname)

}else if (plot.patch){
   cat(" - Retrieving patch-level summary...","\n")
   dummy = load(summary.fullname)
}#end if (! file.exists(summary.fullname))
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Site-specific plots.                                                                 #
#------------------------------------------------------------------------------------------#
if (plot.site){
   #---------------------------------------------------------------------------------------#
   #    Site loop.                                                                         #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(nsites)){
      #----- Retrieve site information. ---------------------------------------------------#
      iata     = sites$iata[p]
      model    = eft[[iata]]
      yeara    = model$yeara
      yearz    = model$yearz
      lon      = model$lon
      lat      = model$lat
      desc     = model$longname
      cat(" - Generating profiles for ",desc,"...","\n")
      #------------------------------------------------------------------------------------#


      #----- Make output path for this place. ---------------------------------------------#
      sitepath = file.path(metroot,model$short)
      if (! file.exists(sitepath)) dir.create(sitepath)
      #------------------------------------------------------------------------------------#



      #----- Choose the patches to plot based on the LAI. ---------------------------------#
      cat("   > LAI by month...","\n")
      if (plot.month){
         lai.mon  = qapply(X = model$lai.pa, INDEX  = model$month, DIM = 2, FUN = mean)
         lai.mon  = apply (X = lai.mon     , MARGIN = c(2,3), FUN = mean)
         loop.mon = sequence(12)
      }else{
         loop.mon = numeric(0)
      }#end if (plot.month)
      #------------------------------------------------------------------------------------#



      #----- Select patches to be plotted. ------------------------------------------------#
      lai.bar      = apply (X=model$lai.pa,MARGIN=3,FUN=mean,na.rm=TRUE)
      use.lai.plot = c( ref = which.closest( x    = sites$lai[p]
                                           , A    = lai.bar
                                           , mask = lai.bar-sites$lai[p]>0
                                           )#end which.closest
                      , max = which.max(lai.bar)
                      , mid = which.closest(x=median(lai.bar),A=lai.bar)
                      , min = which.min(lai.bar)
                      )#end c
      use.lai.key  = c("lairef","laimax","laimed","laimin")
      use.lai.desc = c("Default","Greenest","Medium","Brownest")
      if (plot.refonly){
         use.lai.plot = use.lai.plot[1]
         use.lai.key  = use.lai.key [1]
         use.lai.desc = use.lai.desc[1]
      }#end if (plot.refonly)
      #------------------------------------------------------------------------------------#



      #----- Make output path for this place. ---------------------------------------------#
      profilepath  = file.path(sitepath,"profile" )
      dielpath     = file.path(sitepath,"diel"    )
      reldielpath  = file.path(sitepath,"rel_diel")
      laipath      = file.path(sitepath,"lai"     )
      if (! file.exists(profilepath)) dir.create(profilepath)
      if (! file.exists(dielpath   )) dir.create(dielpath   )
      if (! file.exists(reldielpath)) dir.create(reldielpath)
      if (! file.exists(laipath    )) dir.create(laipath    )
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Loop over the LAI profiles we will plot.                                       #
      #------------------------------------------------------------------------------------#
      for (l in seq_along(use.lai.plot)){
         pnow      = use.lai.plot[l]
         lai.key   = use.lai.key [l]
         lai.desc  = use.lai.desc[l]
         pp        = model$ppaco == pnow
         ptb       = c(which.min(which(pp)),which.max(which(pp)))

         cat("   * ",lai.desc," patch","\n")


         #----- Loop over months and hours. -----------------------------------------------#
         for (m in loop.mon){
            mm   = sprintf("%2.2i",m)
            ms   = model$month %in% m
            cat("     > ",month.name[m],"\n")
            cat("       ~ Average profile by hour of the day...","\n")


            for (h in prof.hours){
               hh      = sprintf("%2.2i",h)
               hs      = model$hour %in% h
               zen.bar = mean(model$zen[ms,hs])



               #------ Create path for this hour. -----------------------------------------#
               profhrpath  = file.path(profilepath,paste(hh,"utc",sep=""))
               if (! file.exists(profhrpath)) dir.create(profhrpath)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Load LAI/WAI/TAI profiles.                                            #
               #---------------------------------------------------------------------------#
               lai     = apply(X=model$lai[,ms,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               wai     = apply(X=model$wai[,ms,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               tai     = apply(X=model$tai[,ms,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               cum.lai = t(apply(X=lai,MARGIN=1,FUN=cumsum))
               cum.wai = t(apply(X=wai,MARGIN=1,FUN=cumsum))
               cum.tai = t(apply(X=tai,MARGIN=1,FUN=cumsum))
               #---------------------------------------------------------------------------#


               #----- Fix the limits for the Y axis. --------------------------------------#
               lai.lim = range(-cum.lai)
               lai.at   = pretty(lai.lim)
               lai.lab  = -lai.at
               #---------------------------------------------------------------------------#



               #----- Band loop. ----------------------------------------------------------#
               for (b in sequence(nband)){
                  bnd      = band$key [b]
                  bnd.desc = band$desc[b]

                  mult     = if (bnd %in% c("par")){Watts.2.Ein * 1.e6}else{1.0}
                  prof     = model[[bnd]]

                  #----- Find the averages for the hour and month. ------------------------#
                  down.beam = apply( X      = prof$down.beam[,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  down.diff = apply( X      = prof$down.diff[,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  down.tot  = apply( X      = prof$down.tot [,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  up.tot    = apply( X      = prof$up.tot   [,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  abs.lyr   = apply( X      = prof$abs.lyr  [,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  abs.tot   = apply( X      = prof$abs.tot  [,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  abs.leaf  = apply( X      = prof$abs.leaf [,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  abs.wood  = apply( X      = prof$abs.wood [,ms,hs,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,4)
                                   , FUN    = mean
                                   )#end apply
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Shift the results a bit so we can see when two curves are on top   #
                  # of each other.                                                         #
                  #------------------------------------------------------------------------#
                  if ( ! bnd %in% "tir"){
                     down.beam = apply(X = down.beam,MARGIN=2,FUN="*",simul$mult.beam.down)
                     down.diff = apply(X = down.diff,MARGIN=2,FUN="*",simul$mult.diff.down)
                     up.tot    = apply(X = up.tot   ,MARGIN=2,FUN="*",simul$mult.up       )
                  }#end if
                  #------------------------------------------------------------------------#





                  #------------------------------------------------------------------------#
                  #     Skip plot if it is solar radiation and it is nighttime.            #
                  #------------------------------------------------------------------------#
                  if (max(down.tot) > 1.){


                     #----- Retrieve the data. --------------------------------------------#
                     if (! bnd %in% "tir"){
                        down.lim = pretty.xylim(u=c(down.beam,down.diff))
                     }else{
                        down.lim = pretty.xylim(u=down.diff)
                     }#end if
                     up.lim   = pretty.xylim(u=up.tot )
                     abs.lim  = pretty.xylim(u=c(0,abs.tot,abs.leaf))
                     #---------------------------------------------------------------------#


                     for (o in sequence(nout)){
                        #----- Open file. -------------------------------------------------#
                        fichier = file.path(profhrpath
                                           ,paste0("profile-",iata,"-",bnd,"-",imetrad,"-"
                                                  ,mm,"-",hh,"utc-",lai.key,".",outform[o]))
                        if (bnd %in% "tir"){nsize=wsize}else{nsize=size}
                        if       (outform[o] == "x11"){
                          X11(width=nsize$width,height=nsize$height,pointsize=ptsz)
                        }else if (outform[o] == "png"){
                          png(filename=fichier,width=nsize$width*depth
                             ,height=nsize$height*depth,pointsize=ptsz,res=depth)
                        }else if (outform[o] == "eps"){
                          postscript(file=fichier,width=nsize$width,height=nsize$height
                                    ,pointsize=ptsz,paper=nsize$paper)
                        }else if (outform[o] == "pdf"){
                          pdf(file=fichier,onefile=FALSE,width=nsize$width
                             ,height=nsize$height,pointsize=ptsz,paper=nsize$paper)
                        }#end if
                        #------------------------------------------------------------------#


                        #----- Set up the plotting window. --------------------------------#
                        par(par.user)
                        par(oma=c(0,1,4,0))
                        if (bnd %in% "tir"){
                           layout(mat=rbind(c(2,3,4),c(1,1,1)),heights=c(6,1))
                        }else{
                           layout(mat=rbind(c(2,4),c(3,5),c(1,1)),heights=c(3,3,1))
                        }#end if
                        #------------------------------------------------------------------#




                        #----- First plot: the legend. ------------------------------------#
                        par(mar=c(0.1,4.1,0.1,2.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1))
                        legend( x      = "bottom"
                              , inset  = 0.0
                              , legend = simul$desc
                              , col    = simul$colour
                              , lty    = "solid"
                              , lwd    = 2.5
                              , cex    = 0.9 * cex.ptsz
                              , ncol   = 2
                              , xpd    = TRUE
                              )#end legend
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Direct radiation is not plotted for TIR.                    #
                        #------------------------------------------------------------------#
                        if (! bnd %in% "tir"){
                           #------ Second plot: Downward direct. --------------------------#
                           par(mar=c(3.1,3.1,3.1,1.1))
                           plot.new()
                           plot.window(xlim=down.lim,ylim=lai.lim)
                           axis(side=1,cex.axis=0.8,padj=-0.8)
                           axis(side=2,at=lai.at,labels=lai.lab,las=1)
                           grid(col=grid.colour,lwd=1.0,lty="solid")
                           box()
                           title(main="Downward - Direct",font.main=1,line=1)
                           for (s in sequence(nsimul)){
                              lines(x=down.beam[s,],y=-cum.lai[s,]
                                   ,col=simul$colour[s],lwd=2.5,lty="solid")
                           }#end for
                           #---------------------------------------------------------------#
                        }#end if
                        #------------------------------------------------------------------#




                        #------ Third plot: Downward diffuse. -----------------------------#
                        par(mar=c(4.1,3.1,2.1,1.1))
                        plot.new()
                        plot.window(xlim=down.lim,ylim=lai.lim)
                        axis(side=1,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=lai.at,labels=lai.lab,las=1)
                        grid(col=grid.colour,lwd=1.0,lty="solid")
                        box()
                        title(main="Downward - Diffuse",font.main=1,line=1)
                        for (s in sequence(nsimul)){
                           lines(x=down.diff[s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Fourth plot: Upward radiation.  Margins depend on the band.  #
                        #------------------------------------------------------------------#
                        if (bnd %in% "tir"){
                           par(mar=c(4.1,3.1,2.1,1.1))
                        }else{
                           par(mar=c(3.1,2.1,3.1,2.1))
                        }#end if
                        plot.new()
                        plot.window(xlim=up.lim,ylim=lai.lim)
                        axis(side=1,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=lai.at,labels=lai.lab,las=1)
                        grid(col=grid.colour,lwd=1.0,lty="solid")
                        box()
                        title(main="Upward",font.main=1,line=1)
                        for (s in sequence(nsimul)){
                           lines(x=up.tot[s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#



                        #------ Fifth plot: Absorbed radiation. ---------------------------#
                        par(mar=c(4.1,2.1,2.1,2.1))
                        plot.new()
                        plot.window(xlim=abs.lim,ylim=lai.lim)
                        axis(side=1,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=lai.at,labels=lai.lab,las=1)
                        grid(col=grid.colour,lwd=1.0,lty="solid")
                        box()
                        title(main="Cumulative absorption (dashed = leaves)"
                             ,font.main=1,line=1)
                        for (s in sequence(nsimul)){
                           lines(x=abs.leaf[s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.0,lty="dotdash")
                           lines(x=abs.tot [s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Plot the general annotation.                                 #
                        #------------------------------------------------------------------#
                        lex     = desc.unit( desc = paste(band$desc[b],"Irradiance")
                                           , unit = untab[[band$unit[b]]]
                                           )#end desc.unit
                        ley     = desc.unit(desc="Cumulative LAI",unit=untab$m2lom2)
                        whenlab = paste(month.name[m],"-",hh,"UTC")
                        letitre = paste(desc," - Canopy profiles","  -  ",whenlab,sep="")
                        lelai   = desc.unit( desc = paste("Average LAI"
                                                         ,sprintf("%6.2f",lai.mon[m,pnow])
                                                         )#end paste
                                           , unit = untab$m2lom2
                                           , bracket = FALSE
                                           )#end desc.unit
                        lezen   = desc.unit( desc    = paste("Sun zenith angle"
                                                            ,sprintf("%5.2f",zen.bar)
                                                            )#end paste
                                           , unit    = untab$deg
                                           , bracket = FALSE
                                           )#end desc.unit
                        lesub   = parse(text = paste(lelai,lezen,sep=" - "))
                        gtitle( main      = letitre
                              , sub       = lesub
                              , xlab      = lex
                              , ylab      = ley
                              , line.xlab = 3.2
                              , line.ylab = 2.6
                              , line.sub  = 4.2
                              , off.sub   = 2/21
                              , cex.main  = 1.0*cex.ptsz
                              , cex.xlab  = 0.9*cex.ptsz
                              , cex.ylab  = 0.9*cex.ptsz
                              , cex.sub   = 0.7*cex.ptsz
                              )#end gtitle
                        #------------------------------------------------------------------#



                        #----- Close the device. ------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        clean.tmp()
                        #------------------------------------------------------------------#
                     }#end for (o in sequence(nout))
                     #---------------------------------------------------------------------#
                  }#end if (max(tot.down) > 1.)
                  #------------------------------------------------------------------------#
               }#end for (b in sequence(nband))
               #---------------------------------------------------------------------------#
            }#end for (h in prof.hours)
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #      Create diel plots for radiation profile at the bottom and top of the    #
            # canopy.                                                                      #
            #------------------------------------------------------------------------------#
            cat("       ~ Diel plots above and beneath canopy.","\n")

            mn = model$month == m

            #----- Initialise arrays that will hold the mean diel. ------------------------#
            empty     = array( data     = NA
                             , dim      = c(nband,nsimul,nhours,2)
                             , dimnames = list(band$key,simul$suff
                                              ,my.hours,c("top","bottom"))
                             )#end array
            flux      = list( down.beam = empty
                            , down.diff = empty
                            , down.tot  = empty
                            , up.tot    = empty
                            )#end list
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Find the mean diurnal flag for each hour.  We will skip the nighttime   #
            # for relative components, to avoid weird results.                             #
            #------------------------------------------------------------------------------#
            sel.day = apply(X=model$diel[mn,,drop=FALSE],MARGIN=2,FUN=commonest)
            sel.day = sel.day %==% 2
            #------------------------------------------------------------------------------#



            #----- Band loop. -------------------------------------------------------------#
            for (b in sequence(nband)){
               bnd      = band$key [b]
               bnd.desc = band$desc[b]
               prof     = model[[bnd]]


               #----- Find the averages for the hour and month. ---------------------------#
               flux$down.beam[b,,,] = apply( X      = prof$down.beam[,mn,,ptb,drop=FALSE]
                                           , MARGIN = c(1,3,4)
                                           , FUN    = mean
                                           )#end apply
               flux$down.diff[b,,,] = apply( X      = prof$down.diff[,mn,,ptb,drop=FALSE]
                                           , MARGIN = c(1,3,4)
                                           , FUN    = mean
                                           )#end apply
               flux$down.tot [b,,,] = apply( X      = prof$down.tot [,mn,,ptb,drop=FALSE]
                                           , MARGIN = c(1,3,4)
                                           , FUN    = mean
                                           )#end apply
               flux$up.tot   [b,,,] = apply( X      = prof$up.tot   [,mn,,ptb,drop=FALSE]
                                           , MARGIN = c(1,3,4)
                                           , FUN    = mean
                                           )#end apply
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Define range for x and y.                                                #
            #------------------------------------------------------------------------------#
            hour.limit = range(my.hours)
            hour.at    = prof.hours
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Loop over the components.  Since shortwave and longwave are very        #
            # different, we only plot shortwave components, and we only use the hours that #
            # are "sufficiently diurnal".                                                  #
            #------------------------------------------------------------------------------#
            nflux = length(flux)
            for (f in sequence(nflux)){
               #----- Get this flux component. --------------------------------------------#
               flux.key   = names(flux)[f]
               flux.desc  = c("Downward (direct)","Downward (diffuse)"
                             ,"Downward (total)" ,"Upward")[f]
               if (flux.key %in% c("down.beam","down.diff")){
                  bwhich     = which(! band$key %in% "tir")
                  nbwhich    = length(bwhich)
               }else{
                  bwhich     = sequence(nband)
                  nbwhich    = nband
               }#end if
               this       = flux[[f]]
               this.limit = range(c(this[bwhich,,,]))
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Loop over formats.                                                   #
               #---------------------------------------------------------------------------#
               for (o in sequence(nout)){
                  #----- Open file. -------------------------------------------------------#
                  fichier = file.path(dielpath
                                     ,paste0("diel-",iata,"-",flux.key,"-",imetrad,"-",mm
                                            ,"-",lai.key,".",outform[o]))
                  if       (outform[o] == "x11"){
                    X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if (outform[o] == "png"){
                    png(filename=fichier,width=size$width*depth
                       ,height=size$height*depth,pointsize=ptsz,res=depth)
                  }else if (outform[o] == "eps"){
                    postscript(file=fichier,width=size$width,height=size$height
                              ,pointsize=ptsz,paper=size$paper)
                  }else if (outform[o] == "pdf"){
                    pdf(file=fichier,onefile=FALSE,width=size$width
                       ,height=size$height,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Set up the plotting window. --------------------------------------#
                  par(par.user)
                  par(oma=c(0,1,4,0))
                  layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                        , heights = c(rep(6/lo.simul$nrow,lo.simul$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      First plot: the legend.  No reason to plot thermal IR radiation   #
                  # for direct component.                                                  #
                  #------------------------------------------------------------------------#
                  if (flux.key %in% c("down.beam","down.diff")){
                     #------ Select the bands that should be used. ------------------------#
                     par(mar=c(0.1,4.1,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = 0.0
                           , legend = c(band$desc[bwhich],"","Top","Bottom")
                           , fill   = c(band$col [bwhich],rep("transparent",3))
                           , border = c(band$col [bwhich],rep("transparent",3))
                           , col    = c(rep("transparent",nband),rep(foreground,2))
                           , lty    = c(rep("solid",nband),"solid","dashed")
                           , lwd    = 2.5
                           , cex    = 0.9 * cex.ptsz
                           , ncol   = 3
                           , xpd    = TRUE
                           )#end legend
                     #---------------------------------------------------------------------#
                  }else{
                     #------ All bands shall be used. -------------------------------------#
                     par(mar=c(0.1,4.1,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = 0.0
                           , legend = c(band$desc,"Top","Bottom")
                           , fill   = c(band$col ,rep("transparent",2))
                           , border = c(band$col ,rep("transparent",2))
                           , col    = c(rep("transparent",nband),rep(foreground,2))
                           , lty    = c(rep("solid",nband),"solid","dashed")
                           , lwd    = 2.5
                           , cex    = 0.9 * cex.ptsz
                           , ncol   = 3
                           , xpd    = TRUE
                           )#end legend
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Loop over the simulations, and the bands.                         #
                  #------------------------------------------------------------------------#
                  for (s in sequence(nsimul)){
                     #------ Open plotting device. ----------------------------------------#
                     par(mar=lo.simul$mar0)
                     plot.new()
                     plot.window(xlim = hour.limit,ylim = this.limit)
                     axis(side=1,at=hour.at)
                     axis(side=2,las=1)
                     abline(v=hour.at,h=axTicks(2),col=grid.colour,lty="dotted")
                     #---------------------------------------------------------------------#



                     #----- Loop over bands to be plotted. --------------------------------#
                     for (b in bwhich){
                         lines(x=my.hours,y=this[b,s,,1],col=band$col[b],type="o",pch=16
                              ,lty="solid" ,lwd=2.0)
                         lines(x=my.hours,y=this[b,s,,2],col=band$col[b],type="o",pch=16
                              ,lty="dashed",lwd=2.0)
                     }#end for (b in bwhich)
                     #---------------------------------------------------------------------#


                     #----- Plot the final information. -----------------------------------#
                     box()
                     title(main=simul$desc[s],font.main=1,line=1)
                     #---------------------------------------------------------------------#
                  }#end for (s in sequence(nsimul))
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot the general annotation.                                       #
                  #------------------------------------------------------------------------#
                  lex     = desc.unit(desc="Hour"      ,unit=untab$utc )
                  ley     = desc.unit(desc="Irradiance",unit=untab$wom2)
                  lesub   = desc.unit( desc    = paste("Average LAI"
                                                      ,sprintf("%6.2f",lai.mon[m,pnow])
                                                      )#end paste
                                     , unit    = untab$m2lom2
                                     , bracket = FALSE
                                     )#end desc.unit
                  letitre = paste(desc," - Mean diel of ",flux.desc,"  -  ",month.name[m]
                                 ,sep="")
                  gtitle( main      = letitre
                        , sub       = lesub
                        , xlab      = lex
                        , ylab      = ley
                        , line.xlab = 3.2
                        , line.ylab = 2.6
                        , line.sub  = 4.2
                        , off.sub   = 2/21
                        , cex.main  = 1.0*cex.ptsz
                        , cex.xlab  = 0.9*cex.ptsz
                        , cex.ylab  = 0.9*cex.ptsz
                        , cex.sub   = 0.7*cex.ptsz
                        )#end gtitle
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in sequence(nout))
               #---------------------------------------------------------------------------#

            }#end for (f in sequence(nflux))
            #------------------------------------------------------------------------------#
         }#end for (m in sequence(12))
         #---------------------------------------------------------------------------------#
      }#end for (l in seq_along(model$lai.plot)
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
#------------------------------------------------------------------------------------------#
}#end if (plot.site)
#==========================================================================================#
#==========================================================================================#
